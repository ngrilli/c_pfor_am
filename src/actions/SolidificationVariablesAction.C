// Nicolò Grilli
// Università di Bristol
// 28 Agosto 2023

#include "SolidificationVariablesAction.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "Factory.h"
#include "FEProblem.h"
#include "NonlinearSystemBase.h"

#include "libmesh/string_to_enum.h"

registerMooseAction("PhaseFieldApp", SolidificationVariablesAction, "add_variable");
registerMooseAction("PhaseFieldApp", SolidificationVariablesAction, "copy_nodal_vars");
registerMooseAction("PhaseFieldApp", SolidificationVariablesAction, "check_copy_nodal_vars");

InputParameters
SolidificationVariablesAction::validParams()
{
  InputParameters params = PolycrystalVariablesAction::validParams();
  params.addClassDescription("Set up order parameter variables for a polycrystal simulation, "
                             "including the zeta variable for solidification.");
  return params;
}

SolidificationVariablesAction::SolidificationVariablesAction(const InputParameters & params)
  : PolycrystalVariablesAction(params),
    _op_num_solidification(getParam<unsigned int>("op_num")),
    _var_name_base_solidification(getParam<std::string>("var_name_base"))
{
}

void
SolidificationVariablesAction::act()
{
  // take initial values from file?
  bool initial_from_file = getParam<bool>("initial_from_file");

  // Loop through the number of order parameters
  for (unsigned int op = 0; op < _op_num_solidification; op++)
  {
    // Create variable names
    std::string var_name = _var_name_base_solidification + Moose::stringify(op);

    // Add the variable
    if (_current_task == "add_variable")
    {
      auto fe_type = AddVariableAction::feType(_pars);
      auto type = AddVariableAction::variableType(fe_type);
      auto var_params = _factory.getValidParams(type);

      var_params.applySpecificParameters(_pars, {"order", "family", "block"});
      var_params.set<std::vector<Real>>("scaling") = {_pars.get<Real>("scaling")};
      _problem->addVariable(type, var_name, var_params);
    }

    // Setup initial from file if requested
    if (initial_from_file)
    {
      if (_current_task == "check_copy_nodal_vars")
        _app.setExodusFileRestart(true);

      if (_current_task == "copy_nodal_vars")
      {
        auto * system = &_problem->getNonlinearSystemBase();
        system->addVariableToCopy(var_name, var_name, "LATEST");
      }
    }
  }
  
  // Create zeta solidification variable
  std::string zeta_name = "zeta";  

  if (_current_task == "add_variable") {
	  
    auto fe_type = AddVariableAction::feType(_pars);
    auto type = AddVariableAction::variableType(fe_type);
    auto var_params = _factory.getValidParams(type);

    var_params.applySpecificParameters(_pars, {"order", "family", "block"});
    var_params.set<std::vector<Real>>("scaling") = {_pars.get<Real>("scaling")};
    _problem->addVariable(type, zeta_name, var_params);    
	  
  }

  if (initial_from_file) {
	  
    if (_current_task == "check_copy_nodal_vars")
      _app.setExodusFileRestart(true);

    if (_current_task == "copy_nodal_vars")
    {
      auto * system = &_problem->getNonlinearSystemBase();
      system->addVariableToCopy(zeta_name, zeta_name, "LATEST");
    }
	  
  }
}
