// Nicolò Grilli
// University of Bristol
// 10 Marzo 2022

#include "CrystalPlasticityStateVarRateBackstress.h"
#include <cmath>

#include <fstream>

registerMooseObject("TensorMechanicsApp", CrystalPlasticityStateVarRateBackstress);

InputParameters
CrystalPlasticityStateVarRateBackstress::validParams()
{
  InputParameters params = CrystalPlasticityStateVarRateComponent::validParams();
  params.addClassDescription("Evolution rate of the backstress based on Armstrong-Frederick. ");  
  params.addParam<std::string>(
      "uo_slip_rate_name",
      "Name of slip rate property: Same as slip rate user object specified in input file.");
  params.addParam<std::string>("uo_backstress_var_name",
                               "Name of backstress variable property: Same as "
                               "backstress variable user object specified in input file. ");
  // Reading the two backstress parameters from file is not implemented
  params.addParam<std::vector<Real>>("bprops","Two parameters required to calculate backstress");
  return params;
}

CrystalPlasticityStateVarRateBackstress::CrystalPlasticityStateVarRateBackstress(const InputParameters & parameters)
  : CrystalPlasticityStateVarRateComponent(parameters),
    _mat_prop_slip_rate(
      getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_slip_rate_name"))),
    _mat_prop_backstress(
      getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_backstress_var_name"))),
    _bprops(getParam<std::vector<Real>>("bprops"))
{
}

bool
CrystalPlasticityStateVarRateBackstress::calcStateVariableEvolutionRateComponent(
    unsigned int qp, std::vector<Real> & val) const
{
  val.assign(_variable_size, 0.0);
  
  // Backstress parameters of the Armstrong-Frederick law
  Real ha = _bprops[0];
  Real hd = _bprops[1];

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    val[i] = ha * _mat_prop_slip_rate[qp][i] 
	       - hd * _mat_prop_backstress[qp][i] * std::abs(_mat_prop_slip_rate[qp][i]);
  }

  return true;
}
