// Nicolò Grilli
// Università di Bristol
// 28 Febbraio 2024

#include "DislocationVariablesAction.h"
#include "FEProblem.h"
#include "Factory.h"

registerMooseAction("c_pfor_amApp", DislocationVariablesAction, "add_aux_variable");
registerMooseAction("c_pfor_amApp", DislocationVariablesAction, "add_aux_kernel");

InputParameters
DislocationVariablesAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Set up auxvariables and auxkernels to visualize dislocation densities");
  params.addParam<std::string>("base_name", "Auxiliary variables base name");
  params.addRequiredParam<unsigned int>(
      "number_slip_systems",
      "The total number of possible active slip systems for the crystalline material");
  params.addParam<bool>("rho_ssd", false, "Add rho_ssd variables");
  params.addParam<bool>("rho_gnd_edge", false, "Add rho_gnd_edge variables");
  params.addParam<bool>("rho_gnd_screw", false, "Add rho_gnd_screw variables");
  params.addParam<bool>("slip_resistance", false, "Add slip_resistance variables");
  return params;
}

DislocationVariablesAction::DislocationVariablesAction(const InputParameters & params)
  : Action(params),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _number_slip_systems(getParam<unsigned int>("number_slip_systems")),
    _rho_ssd(getParam<bool>("rho_ssd")),
    _rho_gnd_edge(getParam<bool>("rho_gnd_edge")),
    _rho_gnd_screw(getParam<bool>("rho_gnd_screw")),
    _slip_resistance(getParam<bool>("slip_resistance"))
{
}

void
DislocationVariablesAction::act()
{
  if (_current_task == "add_aux_variable")
  {
    // Add SSD dislocation density auxvariables
    if (_rho_ssd) {
      for (const auto i : make_range(_number_slip_systems)) {
      
        std::string var_name = _base_name + "rho_ssd_" + Moose::stringify(i);
      
        auto var_params = _factory.getValidParams("MooseVariableConstMonomial");
        var_params.set<MooseEnum>("family") = "MONOMIAL";
        var_params.set<MooseEnum>("order") = "CONSTANT";
      
        _problem->addAuxVariable("MooseVariable", var_name, var_params);
      }
    }
    // Add edge GND dislocation density auxvariables
    if (_rho_gnd_edge) {
      for (const auto i : make_range(_number_slip_systems)) {
      
        std::string var_name = _base_name + "rho_gnd_edge_" + Moose::stringify(i);
      
        auto var_params = _factory.getValidParams("MooseVariableConstMonomial");
        var_params.set<MooseEnum>("family") = "MONOMIAL";
        var_params.set<MooseEnum>("order") = "CONSTANT";
      
        _problem->addAuxVariable("MooseVariable", var_name, var_params);
      }
    }
    // Add screw GND dislocation density auxvariables
    if (_rho_gnd_screw) {
      for (const auto i : make_range(_number_slip_systems)) {
      
        std::string var_name = _base_name + "rho_gnd_screw_" + Moose::stringify(i);
      
        auto var_params = _factory.getValidParams("MooseVariableConstMonomial");
        var_params.set<MooseEnum>("family") = "MONOMIAL";
        var_params.set<MooseEnum>("order") = "CONSTANT";
      
        _problem->addAuxVariable("MooseVariable", var_name, var_params);
      }
    }
    // Add slip resistance auxvariables
    if (_slip_resistance) {
      for (const auto i : make_range(_number_slip_systems)) {
      
        std::string var_name = _base_name + "slip_resistance_" + Moose::stringify(i);
      
        auto var_params = _factory.getValidParams("MooseVariableConstMonomial");
        var_params.set<MooseEnum>("family") = "MONOMIAL";
        var_params.set<MooseEnum>("order") = "CONSTANT";
      
        _problem->addAuxVariable("MooseVariable", var_name, var_params);
      }
    }
  } 
  else if (_current_task == "add_aux_kernel")
  {
	// Add SSD dislocation density auxkernels
	if (_rho_ssd) {
      for (const auto i : make_range(_number_slip_systems)) {

        std::string kernel_name = _base_name + "rho_ssd_" + Moose::stringify(i);
      
        auto kernel_params = _factory.getValidParams("MaterialStdVectorAux");
        kernel_params.set<AuxVariableName>("variable") = kernel_name;
        kernel_params.set<MaterialPropertyName>("property") = "rho_ssd";
        kernel_params.set<unsigned int>("index") = i;
        kernel_params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
        kernel_params.applyParameters(parameters());

        _problem->addAuxKernel("MaterialStdVectorAux", kernel_name, kernel_params);
      }
    }
	// Add edge GND dislocation density auxkernels
	if (_rho_gnd_edge) {
      for (const auto i : make_range(_number_slip_systems)) {

        std::string kernel_name = _base_name + "rho_gnd_edge_" + Moose::stringify(i);
      
        auto kernel_params = _factory.getValidParams("MaterialStdVectorAux");
        kernel_params.set<AuxVariableName>("variable") = kernel_name;
        kernel_params.set<MaterialPropertyName>("property") = "rho_gnd_edge";
        kernel_params.set<unsigned int>("index") = i;
        kernel_params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
        kernel_params.applyParameters(parameters());

        _problem->addAuxKernel("MaterialStdVectorAux", kernel_name, kernel_params);
      }
    }  
	// Add screw GND dislocation density auxkernels
	if (_rho_gnd_screw) {
      for (const auto i : make_range(_number_slip_systems)) {

        std::string kernel_name = _base_name + "rho_gnd_screw_" + Moose::stringify(i);
      
        auto kernel_params = _factory.getValidParams("MaterialStdVectorAux");
        kernel_params.set<AuxVariableName>("variable") = kernel_name;
        kernel_params.set<MaterialPropertyName>("property") = "rho_gnd_screw";
        kernel_params.set<unsigned int>("index") = i;
        kernel_params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
        kernel_params.applyParameters(parameters());

        _problem->addAuxKernel("MaterialStdVectorAux", kernel_name, kernel_params);
      }
    }      
	// Add slip resistance auxkernels
	if (_slip_resistance) {
      for (const auto i : make_range(_number_slip_systems)) {

        std::string kernel_name = _base_name + "slip_resistance_" + Moose::stringify(i);
      
        auto kernel_params = _factory.getValidParams("MaterialStdVectorAux");
        kernel_params.set<AuxVariableName>("variable") = kernel_name;
        kernel_params.set<MaterialPropertyName>("property") = "slip_resistance";
        kernel_params.set<unsigned int>("index") = i;
        kernel_params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
        kernel_params.applyParameters(parameters());

        _problem->addAuxKernel("MaterialStdVectorAux", kernel_name, kernel_params);
      }
    }
  }
}
