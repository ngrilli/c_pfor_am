// Nicolò Grilli
// Università di Bristol
// 10 Dicembre 2023

#include "LiquidSolidKernelAction.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "Factory.h"
#include "FEProblem.h"
#include "NonlinearSystemBase.h"

#include "libmesh/string_to_enum.h"

registerMooseAction("c_pfor_amApp", LiquidSolidKernelAction, "add_variable");
registerMooseAction("c_pfor_amApp", LiquidSolidKernelAction, "add_kernel");

InputParameters
LiquidSolidKernelAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Set up Reaction, BodyForce, CoupledTanh, ..., Diffusion kernels "
      "for the zeta variable, which is 0 in the liquid phase and 1 in the solid phase. ");
      
  // Get MooseEnums for the possible order/family options for the zeta variable
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  params.addParam<MooseEnum>("family",
                             families,
                             "Specifies the family of FE "
                             "shape function to use for the zeta variable");
  params.addParam<MooseEnum>("order",
                             orders,
                             "Specifies the order of the FE "
                             "shape function to use for the zeta variable");
  
  return params;
}

LiquidSolidKernelAction::LiquidSolidKernelAction(const InputParameters & params)
  : Action(params)  
{
}

void
LiquidSolidKernelAction::act()
{
  std::string zeta_var_name = "zeta";
  
  // Add the zeta variable
  if (_current_task == "add_variable")
  { 
    auto fe_type = AddVariableAction::feType(_pars);
    auto type = AddVariableAction::variableType(fe_type);
    auto var_params = _factory.getValidParams(type);

    var_params.applySpecificParameters(_pars, {"order", "family", "block"});
    var_params.set<std::vector<Real>>("scaling") = {_pars.get<Real>("scaling")};
    _problem->addVariable(type, zeta_var_name, var_params);
  }	
  
  if (_current_task == "add_kernel") {
    // Add the reaction kernel
    InputParameters params = _factory.getValidParams("Reaction");
    params.set<NonlinearVariableName>("variable") = zeta_var_name;
    // define rate which is the prefactor
    params.applyParameters(parameters());

    std::string kernel_name = "Reaction_" + zeta_var_name;
    _problem->addKernel("Reaction", kernel_name, params);
  }
}
