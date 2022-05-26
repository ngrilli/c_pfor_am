// Nicolò Grilli
// University of Bristol
// 17 Marzo 2022

#include "CrystalPlasticitySlipRateGSSBackstress.h"

#include <fstream>

registerMooseObject("TensorMechanicsApp", CrystalPlasticitySlipRateGSSBackstress);

InputParameters
CrystalPlasticitySlipRateGSSBackstress::validParams()
{
  InputParameters params = CrystalPlasticitySlipRateGSS::validParams();
  params.addParam<std::string>("backstress_var_name",
                               "Name of backstress property: Same as "
                               "state variable user object for backstress specified " 
							   "in the input file. ");
  params.addClassDescription("Phenomenological constitutive model slip rate class "
                             "with backstress. Override the "
                             "virtual functions in your class");
  return params;
}

CrystalPlasticitySlipRateGSSBackstress::CrystalPlasticitySlipRateGSSBackstress(const InputParameters & parameters)
  : CrystalPlasticitySlipRateGSS(parameters),
    _mat_prop_backstress(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("backstress_var_name")))

{
  if (_slip_sys_flow_prop_file_name.length() != 0)
    readFileFlowRateParams();
  else
    getFlowRateParams();
}

bool
CrystalPlasticitySlipRateGSSBackstress::calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);
  }

  // Add backstress to the slip rate law
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    val[i] = _a0(i) * std::pow(std::abs((tau(i) - _mat_prop_backstress[qp][i]) 
	       / _mat_prop_state_var[qp][i]), 1.0 / _xm(i)) *
           std::copysign(1.0, tau(i) - _mat_prop_backstress[qp][i]);

    if (std::abs(val[i] * dt) > _slip_incr_tol)
    {
#ifdef DEBUG
      mooseWarning("Maximum allowable slip increment exceeded ", std::abs(val[i]) * dt);
#endif
      return false;
    }
  }

  return true;
}

bool
CrystalPlasticitySlipRateGSSBackstress::calcSlipRateDerivative(unsigned int qp,
                                                     Real /*dt*/,
                                                     std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    val[i] = _a0(i) / _xm(i) *
           std::pow(std::abs((tau(i) - _mat_prop_backstress[qp][i]) 
		   / _mat_prop_state_var[qp][i]), 1.0 / _xm(i) - 1.0) /
           _mat_prop_state_var[qp][i];
  }
  
  return true;
}
