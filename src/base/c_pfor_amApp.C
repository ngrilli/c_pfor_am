#include "c_pfor_amApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
c_pfor_amApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

c_pfor_amApp::c_pfor_amApp(InputParameters parameters) : MooseApp(parameters)
{
  c_pfor_amApp::registerAll(_factory, _action_factory, _syntax);
}

c_pfor_amApp::~c_pfor_amApp() {}

static void
associateSyntaxInner(Syntax & syntax, ActionFactory & /*action_factory*/)
{
  registerSyntax("PolycrystalSolidificationKernelAction", "Kernels/PolycrystalSolidificationKernel");
  registerSyntax("LiquidSolidKernelAction", "Kernels/LiquidSolidKernel");
  registerSyntax("SolidificationVariablesAction", "Variables/SolidificationVariables");
  registerSyntax("DislocationSlipGradientAction", "AuxVariables/DislocationSlipGradient");
  registerSyntax("DislocationSlipGradientAction", "AuxKernels/DislocationSlipGradient");
  registerSyntax("DislocationVariablesAction", "AuxVariables/DislocationVariables");
  registerSyntax("DislocationVariablesAction", "AuxKernels/DislocationVariables");
}

void
c_pfor_amApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<c_pfor_amApp>(f, af, s);
  Registry::registerObjectsTo(f, {"c_pfor_amApp"});
  Registry::registerActionsTo(af, {"c_pfor_amApp"});

  /* register custom execute flags, action syntax, etc. here */
  associateSyntaxInner(s, af);
}

void
c_pfor_amApp::registerApps()
{
  registerApp(c_pfor_amApp);
}

//void
//c_pfor_amApp::registerObjects(Factory & factory)
//{
//  mooseDeprecated("use registerAll instead of registerObjects");
//  Registry::registerObjectsTo(factory, {"c_pfor_amApp"});
//}

//void
//c_pfor_amApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
//{
//  mooseDeprecated("use registerAll instead of associateSyntax");
//  Registry::registerActionsTo(action_factory, {"c_pfor_amApp"});
//  associateSyntaxInner(syntax, action_factory);
//}

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
