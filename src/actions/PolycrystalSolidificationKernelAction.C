// Nicolò Grilli
// Università di Bristol
// 26 Agosto 2023

#include "PolycrystalSolidificationKernelAction.h"
#include "Factory.h"
#include "Conversion.h"
#include "FEProblem.h"

registerMooseAction("c_pfor_amApp", PolycrystalSolidificationKernelAction, "add_kernel");

InputParameters
PolycrystalSolidificationKernelAction::validParams()
{
  InputParameters params = PolycrystalKernelAction::validParams();
  params.addClassDescription(
      "Set up ACGrGrPoly, GrainSolidification, ACInterface, TimeDerivative, and ACGBPoly kernels");
  params.addParam<VariableName>("zeta", "Phase field representing liquid (0) or solid (1)");
  return params;
}

PolycrystalSolidificationKernelAction::PolycrystalSolidificationKernelAction(const InputParameters & params)
  : PolycrystalKernelAction(params)  
{
}

void
PolycrystalSolidificationKernelAction::act()
{
  // Create zeta solid-liquid variable name
  // The variable zeta must be added to the input file
  // or created by another action
  VariableName zeta = "zeta";
  
  for (unsigned int op = 0; op < _op_num; ++op)
  {
    //
    // Create variable names
    //

    std::string var_name = _var_name_base + Moose::stringify(op);
    std::vector<VariableName> v;
    v.resize(_op_num - 1);

    unsigned int ind = 0;
    for (unsigned int j = 0; j < _op_num; ++j)
      if (j != op)
        v[ind++] = _var_name_base + Moose::stringify(j);

    //
    // Set up ACGrGrPoly kernels
    //

    {
      InputParameters params = _factory.getValidParams("ACGrGrPoly");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<VariableName>>("v") = v;
      params.applyParameters(parameters());

      std::string kernel_name = "ACBulk_" + var_name;
      _problem->addKernel("ACGrGrPoly", kernel_name, params);
    }
    
    //
    // Set up GrainSolidification kernels
    //
    
    if (isParamValid("zeta"))
    {
      InputParameters params = _factory.getValidParams("GrainSolidification");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<VariableName>>("v") = v;
      params.set<VariableName>("zeta") = zeta;
      params.applyParameters(parameters());
      
      std::string kernel_name = "GrainSolidification_" + var_name;
      _problem->addKernel("GrainSolidification", kernel_name, params);
      
	} else {
		
      mooseError("zeta field liquid-solid variable missing in PolycrystalSolidificationAction");
	}

    //
    // Set up ACInterface kernels
    //

    {
      InputParameters params = _factory.getValidParams("ACInterface");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.applyParameters(parameters());

      std::string kernel_name = "ACInt_" + var_name;
      _problem->addKernel("ACInterface", kernel_name, params);
    }

    //
    // Set up TimeDerivative kernels
    //

    {
      InputParameters params = _factory.getValidParams("TimeDerivative");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<bool>("implicit") = true;
      params.applyParameters(parameters());

      std::string kernel_name = "IE_" + var_name;
      _problem->addKernel("TimeDerivative", kernel_name, params);
    }

    //
    // Set up optional ACGBPoly bubble interaction kernels
    //

    if (isParamValid("c"))
    {
      InputParameters params = _factory.getValidParams("ACGBPoly");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<VariableName>>("c") = {getParam<VariableName>("c")};
      params.applyParameters(parameters());

      std::string kernel_name = "ACBubInteraction_" + var_name;
      _problem->addKernel("ACGBPoly", kernel_name, params);
    }
  }
}
