// Nicolò Grilli
// Parsa Esmati
// Università di Bristol
// Fernando Valiente Dies
// ANSTO
// 10 Dicembre 2023

#include "LiquidSolidKernelAction.h"
#include "Conversion.h"
#include "Factory.h"
#include "FEProblem.h"
#include "NonlinearSystemBase.h"

#include "libmesh/string_to_enum.h"

registerMooseAction("c_pfor_amApp", LiquidSolidKernelAction, "add_kernel");

InputParameters
LiquidSolidKernelAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Set up TimeDerivative, Reaction, BodyForce, CoupledTanh, Diffusion, CoupledPhaseGrain kernels "
      "for the zeta variable, which is 0 in the liquid phase and 1 in the solid phase. ");

  // Coupled variables: temperature and phase fields
  params.addParam<VariableName>("temperature", "Temperature. ");
  params.addRequiredParam<unsigned int>(
      "op_num", "specifies the total number of grains (deformed + recrystallized) to create. ");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name "
                                                        "of the phase field variables. ");
                                                      
  // Model parameters
  params.addParam<Real>("sigma_p", 0.385e-12, "Liquid-solid interfacial energy (J/micron^2). ");
  params.addParam<Real>("delta_f_p", 0.25, "The maximum height of the barrier in the free energy density"
                                          " between two minima for liquid-solid system. ");
  params.addParam<Real>("l_p", 9.6, "Liquid-solid interfacial width (micron). ");
  params.addParam<Real>("L_p", 1.0, "Kinetic coefficient for phase transformation. ");
  params.addParam<Real>("theta", 21.0, "Theta coefficient for hyperbolic tangent for liquid-solid. ");
  params.addParam<Real>("T_l", 1723.0, "Liquidus temperature (K). ");
  params.addParam<Real>("gamma_p", 1.0, "Interaction coefficient between zeta and eta variables. ");
  
  return params;
}

LiquidSolidKernelAction::LiquidSolidKernelAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _var_name_base(getParam<std::string>("var_name_base")),
    _sigma_p(getParam<Real>("sigma_p")),
    _delta_f_p(getParam<Real>("delta_f_p")),
    _l_p(getParam<Real>("l_p")),
    _L_p(getParam<Real>("L_p")),
    _theta(getParam<Real>("theta")),
    _T_l(getParam<Real>("T_l")),
    _gamma_p(getParam<Real>("gamma_p"))
{
}

void
LiquidSolidKernelAction::act()
{
  std::string zeta_var_name = "zeta";
  
  // Material constant calculations
  _m_p = (3.0/4.0) * _sigma_p / (_delta_f_p * _l_p);
  _k_p = (3.0/4.0) * _sigma_p * _l_p;
  
  // Create phase field variable names
  std::vector<VariableName> phase_fields;
  phase_fields.resize(_op_num);
  
  unsigned int ind = 0;
  
  for (unsigned int j = 0; j < _op_num; ++j) {
	 
    phase_fields[ind++] = _var_name_base + Moose::stringify(j);
	  
  }

  if (_current_task == "add_kernel") {
    // Add the time derivative kernel, LHS term: \dot{\zeta}
    InputParameters params_TimeDerivative = _factory.getValidParams("TimeDerivative");
    params_TimeDerivative.set<NonlinearVariableName>("variable") = zeta_var_name;
    params_TimeDerivative.set<bool>("implicit") = true;
    params_TimeDerivative.applyParameters(parameters());

    std::string kernel_name_TimeDerivative = "TimeDerivative_" + zeta_var_name;
    _problem->addKernel("TimeDerivative", kernel_name_TimeDerivative, params_TimeDerivative); 
	  
    // Add the reaction kernel, LHS term: + 2 m_p L_p \zeta
    InputParameters params_Reaction = _factory.getValidParams("Reaction");
    params_Reaction.set<NonlinearVariableName>("variable") = zeta_var_name;
    params_Reaction.set<Real>("rate") = 2.0 * _m_p * _L_p;
    params_Reaction.applyParameters(parameters());

    std::string kernel_name_Reaction = "Reaction_" + zeta_var_name;
    _problem->addKernel("Reaction", kernel_name_Reaction, params_Reaction);
    
    // Add the body force kernel, LHS term: -m_p L_p
    InputParameters params_BodyForce = _factory.getValidParams("BodyForce");
    params_BodyForce.set<NonlinearVariableName>("variable") = zeta_var_name;
    params_BodyForce.set<Real>("value") = _m_p * _L_p;
    params_BodyForce.applyParameters(parameters());
    
    std::string kernel_name_BodyForce = "BodyForce_" + zeta_var_name;
    _problem->addKernel("BodyForce", kernel_name_BodyForce, params_BodyForce);
    
    // Add the coupled tanh kernel, LHS term: + m_p L_p tanh(theta (T/T_l - 1))
    InputParameters params_CoupledTanh = _factory.getValidParams("CoupledTanh");
    params_CoupledTanh.set<NonlinearVariableName>("variable") = zeta_var_name;
    params_CoupledTanh.set<std::vector<VariableName>>("v") = {getParam<VariableName>("temperature")};
    params_CoupledTanh.set<Real>("A") = _m_p * _L_p;
    params_CoupledTanh.set<Real>("theta") = _theta;
    params_CoupledTanh.set<Real>("vn") = _T_l;
    params_CoupledTanh.applyParameters(parameters());
    
    std::string kernel_name_CoupledTanh = "CoupledTanh_" + zeta_var_name;
    _problem->addKernel("CoupledTanh", kernel_name_CoupledTanh, params_CoupledTanh);
    
    // Add the diffusion kernel, LHS term: - L_p k_p \nabla^2 \zeta
    InputParameters params_Diffusion = _factory.getValidParams("CoefDiffusion");
    params_Diffusion.set<NonlinearVariableName>("variable") = zeta_var_name;
    params_Diffusion.set<Real>("coef") = _L_p * _k_p;
    params_Diffusion.applyParameters(parameters());
    
    std::string kernel_name_Diffusion = "Diffusion_" + zeta_var_name;
    _problem->addKernel("CoefDiffusion", kernel_name_Diffusion, params_Diffusion);
    
    // Add the coupled phase and grain kernel, LHS term: -2 L_p m_g gamma_p (1 - \zeta) \sum \eta_i^2
    InputParameters params_CoupledPhaseGrain = _factory.getValidParams("CoupledPhaseGrain");
    params_CoupledPhaseGrain.set<NonlinearVariableName>("variable") = zeta_var_name;
    params_CoupledPhaseGrain.set<std::vector<VariableName>>("v") = phase_fields;
    params_CoupledPhaseGrain.set<Real>("L_p") = _L_p;
    params_CoupledPhaseGrain.set<Real>("gamma_p") = _gamma_p;
    params_CoupledPhaseGrain.applyParameters(parameters());
    
    std::string kernel_name_CoupledPhaseGrain = "CoupledPhaseGrain_" + zeta_var_name;
    _problem->addKernel("CoupledPhaseGrain", kernel_name_CoupledPhaseGrain, params_CoupledPhaseGrain);
    
  }
}
