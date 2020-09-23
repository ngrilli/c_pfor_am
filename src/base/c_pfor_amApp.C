#include "c_pfor_amApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
c_pfor_amApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

c_pfor_amApp::c_pfor_amApp(InputParameters parameters) : MooseApp(parameters)
{
  c_pfor_amApp::registerAll(_factory, _action_factory, _syntax);
}

c_pfor_amApp::~c_pfor_amApp() {}

void
c_pfor_amApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"c_pfor_amApp"});
  Registry::registerActionsTo(af, {"c_pfor_amApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
c_pfor_amApp::registerApps()
{
  registerApp(c_pfor_amApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
c_pfor_amApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  c_pfor_amApp::registerAll(f, af, s);
}
extern "C" void
c_pfor_amApp__registerApps()
{
  c_pfor_amApp::registerApps();
}
