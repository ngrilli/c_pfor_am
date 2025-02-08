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
  // Add slip rate auxvariables
  if (_current_task == "add_aux_variable")
  {
    for (const auto i : make_range(_number_slip_systems)) {
      
      std::string var_name = _base_name + "slip_rate_" + Moose::stringify(i);
      
      auto var_params = _factory.getValidParams("MooseVariable");
      var_params.set<MooseEnum>("family") = "MONOMIAL";
      var_params.set<MooseEnum>("order") = "FIRST";
      
      _problem->addAuxVariable("MooseVariable", var_name, var_params);
    }
  } 
  else if (_current_task == "add_aux_kernel")
  {
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
  }
}
