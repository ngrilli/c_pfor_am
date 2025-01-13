// Nicolò Grilli
// Università di Bristol
// 23 Dicembre 2024

#include "DislocationSlipGradientAction.h"
#include "FEProblem.h"
#include "Factory.h"

registerMooseAction("c_pfor_amApp", DislocationSlipGradientAction, "add_aux_variable");

InputParameters
DislocationSlipGradientAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Set up slip rate auxkernels for the dislocation slip gradient model");
  params.addRequiredParam<unsigned int>(
      "number_slip_systems",
      "The total number of possible active slip systems for the crystalline material");
  return params;
}

DislocationSlipGradientAction::DislocationSlipGradientAction(const InputParameters & params)
  : Action(params),
    _number_slip_systems(getParam<unsigned int>("number_slip_systems"))
{
}

void
DislocationSlipGradientAction::act()
{
  // Add slip rate auxvariables
  if (_current_task == "add_aux_variable")
  {
    std::string type = "MooseVariable";

    for (const auto i : make_range(_number_slip_systems)) {
      
      std::string var_name = "slip_rate_" + Moose::stringify(i);
      
      auto var_params = _factory.getValidParams("MooseVariable");
      var_params.set<MooseEnum>("family") = "LAGRANGE";
      var_params.set<MooseEnum>("order") = "FIRST";
      
      _problem->addAuxVariable("MooseVariable", var_name, var_params);
    }
  }
}
