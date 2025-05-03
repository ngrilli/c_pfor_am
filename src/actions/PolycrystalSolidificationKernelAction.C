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
      "Set up ACGrGrPoly, GrainSolidification, ACInterface, TimeDerivative, and ACGBPoly kernels. ");
  params.addParam<VariableName>("zeta", "Phase field representing liquid (0) or solid (1). ");
  params.addParam<Real>("gamma_p", 1.0, "Interaction coefficient between zeta and eta variables. ");
  params.addParam<bool>("GB_anisotropy", false, "Flag to activate grain boundary anisotropy. ");
  params.addParam<Real>("e_anisotropy", 0.0, "Grain boundary energy anisotropy coefficient. ");
  params.addParam<bool>("continuous_anisotropy", false, "Flag to activate continuous model for anisotropy");
  params.addParam<FileName>(
      "Euler_angles_file_name","",
      "Name of the file containing the Euler angles, each row must contain three Euler angles "
      "which correspond to each grain orientation. ");
  return params;
}

PolycrystalSolidificationKernelAction::PolycrystalSolidificationKernelAction(const InputParameters & params)
  : PolycrystalKernelAction(params),
    _gamma_p(getParam<Real>("gamma_p")),
    _GB_anisotropy(this->template getParam<bool>("GB_anisotropy")),
    _e_anisotropy(getParam<Real>("e_anisotropy")),
    _continuous_anisotropy(getParam<bool>("continuous_anisotropy")),
    _Euler_angles_file_name(getParam<FileName>("Euler_angles_file_name"))
{
}

void
PolycrystalSolidificationKernelAction::act()
{
	
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
    // The variable zeta must be added to the input file
    // or created by another action
    //
    
    if (isParamValid("zeta"))
    {
      InputParameters params = _factory.getValidParams("GrainSolidification");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<VariableName>>("v") = v;
      params.set<std::vector<VariableName>>("zeta") = {getParam<VariableName>("zeta")};
      params.set<Real>("gamma_p") = _gamma_p;
      params.applyParameters(parameters());
      
      std::string kernel_name = "GrainSolidification_" + var_name;
      _problem->addKernel("GrainSolidification", kernel_name, params);
      
	} else {
		
      mooseError("zeta field liquid-solid variable missing in PolycrystalSolidificationAction");
	}

    //
    // Set up ACInterface or ACInterfaceAniso kernels
    //

    if (_GB_anisotropy)
    {
      InputParameters params = _factory.getValidParams("ACInterfaceAniso");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<Real>("e_anisotropy") = _e_anisotropy;
      params.set<bool>("continuous_anisotropy") = _continuous_anisotropy;
      params.set<int>("op") = op;
      params.set<int>("op_num") = _op_num;
      
      if (this->isParamValid("Euler_angles_file_name")) {
		  
        params.set<FileName>("Euler_angles_file_name") = _Euler_angles_file_name;
		  
	  }

      params.applyParameters(parameters());

      std::string kernel_name = "ACIntAniso_" + var_name;
      _problem->addKernel("ACInterfaceAniso", kernel_name, params);
      
    } else {

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
