// Nicolò Grilli
// Università di Bristol
// 23 Dicembre 2024

#include "DislocationSlipGradientAction.h"
#include "FEProblem.h"
#include "Factory.h"

registerMooseAction("c_pfor_amApp", DislocationSlipGradientAction, "add_aux_variable");
registerMooseAction("c_pfor_amApp", DislocationSlipGradientAction, "add_aux_kernel");

InputParameters
DislocationSlipGradientAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Set up slip rate auxkernels for the dislocation slip gradient model");
  params.addParam<std::string>("base_name", "Auxiliary variables base name");
  params.addRequiredParam<unsigned int>(
      "number_slip_systems",
      "The total number of possible active slip systems for the crystalline material");
  return params;
}

DislocationSlipGradientAction::DislocationSlipGradientAction(const InputParameters & params)
  : Action(params),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _number_slip_systems(getParam<unsigned int>("number_slip_systems"))
{
}

void
DislocationSlipGradientAction::act()
{
  if (_current_task == "add_aux_variable")
  {
    // Add slip rate auxvariables
    for (const auto i : make_range(_number_slip_systems)) {
      
      std::string var_name = _base_name + "slip_rate_" + Moose::stringify(i);
      
      auto var_params = _factory.getValidParams("MooseVariable");
      var_params.set<MooseEnum>("family") = "MONOMIAL";
      var_params.set<MooseEnum>("order") = "FIRST";
      
      _problem->addAuxVariable("MooseVariable", var_name, var_params);
    }
    // Add slip rate vector auxvariable
    std::string var_name = _base_name + "slip_rate_vector";
    
    auto var_params = _factory.getValidParams("MooseVariable");
    var_params.set<MooseEnum>("family") = "MONOMIAL";
    var_params.set<MooseEnum>("order") = "FIRST";
    var_params.set<unsigned int>("components") = _number_slip_systems;
    
    _problem->addAuxVariable("ArrayMooseVariable", var_name, var_params);
  } 
  else if (_current_task == "add_aux_kernel")
  {
	// Add slip rate auxkernels storing material property slip_increment
    for (const auto i : make_range(_number_slip_systems)) {

      std::string kernel_name = _base_name + "slip_rate_" + Moose::stringify(i);
      
      auto kernel_params = _factory.getValidParams("MaterialStdVectorAux");
      kernel_params.set<AuxVariableName>("variable") = kernel_name;
      kernel_params.set<MaterialPropertyName>("property") = "slip_increment";
      kernel_params.set<unsigned int>("index") = i;
      kernel_params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
      kernel_params.applyParameters(parameters());

      _problem->addAuxKernel("MaterialStdVectorAux", kernel_name, kernel_params);
    }
    // Add auxkernel to build slip rate vector
    std::string kernel_name = _base_name + "slip_rate_vector";
    
    auto kernel_params = _factory.getValidParams("BuildArrayVariableAux");
    kernel_params.set<AuxVariableName>("variable") = kernel_name;
    
    // String with list of slip rate variable names
    std::vector<VariableName> slip_rate_variable_names;
    for (const auto i : make_range(_number_slip_systems)) {
	  slip_rate_variable_names.push_back(_base_name + "slip_rate_" + Moose::stringify(i));
	}
    kernel_params.set<std::vector<VariableName>>("component_variables") = {slip_rate_variable_names};
    
    _problem->addAuxKernel("BuildArrayVariableAux", kernel_name, kernel_params);
  }
}
