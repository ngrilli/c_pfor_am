// Nicolo Grilli
// Daijun Hu 
// National University of Singapore
// 9 Ottobre 2020

#include "FiniteStrainCrystalPlasticityThermal.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

#include <fstream>
#include <cmath>

registerMooseObject("TensorMechanicsApp", FiniteStrainCrystalPlasticityThermal);

InputParameters
FiniteStrainCrystalPlasticityThermal::validParams()
{
  InputParameters params = FiniteStrainCrystalPlasticity::validParams();
  params.addClassDescription("Crystal Plasticity with thermal eigenstrain");
  params.addCoupledVar("temp",293.0,"Temperature");
  params.addParam<Real>("thermal_expansion",0.0,"Thermal expansion coefficient");
  params.addParam<Real>("reference_temperature",293.0,"reference temperature for thermal expansion");
  params.addParam<Real>("dCRSS_dT_A",1.0,"A coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 293.0))");
  params.addParam<Real>("dCRSS_dT_B",0.0,"B coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 293.0))");
  params.addParam<Real>("dCRSS_dT_C",0.0,"C coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 293.0))");
  return params;
}

FiniteStrainCrystalPlasticityThermal::FiniteStrainCrystalPlasticityThermal(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
    _temp(coupledValue("temp")),     
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _dCRSS_dT_A(getParam<Real>("dCRSS_dT_A")),
	_dCRSS_dT_B(getParam<Real>("dCRSS_dT_B")),
	_dCRSS_dT_C(getParam<Real>("dCRSS_dT_C")),
	_gssT(_nss)
{
}

void
FiniteStrainCrystalPlasticityThermal::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;
  Real temp = _temp[_qp];
  Real thermal_expansion = _thermal_expansion; 
  Real reference_temperature = _reference_temperature;


  iden.zero();
  iden.addIa(1.0);

  _fe = _dfgrd_tmp * _fp_prev_inv; // _fp_inv  ==> _fp_prev_inv
  ce = _fe.transpose() * _fe;
  ce_pk2 = ce * _pk2_tmp;
  ce_pk2 = ce_pk2 / _fe.det();

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _nss; ++i)
    _tau(i) = ce_pk2.doubleContraction(_s0[i]);

  getSlipIncrements(); // Calculate dslip,dslipdtau

  if (_err_tol)
    return;

  eqv_slip_incr.zero();
  for (unsigned int i = 0; i < _nss; ++i)
    eqv_slip_incr += _s0[i] * _slip_incr(i);

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;

  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;
  RankTwoTensor thermal_eigenstrain;
  thermal_eigenstrain = (1.0 / 2.0)
                      * (std::exp((2.0/3.0) * thermal_expansion * (temp - reference_temperature)) - 1.0)
                      * iden;
  pk2_new = _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);
  
  resid = _pk2_tmp - pk2_new;
}

// Calculate slip increment,dslipdtau
// Critical resolved shear stress decreases exponentially with temperature
// A + B exp(- C * (T - 293.0))
void
FiniteStrainCrystalPlasticityThermal::getSlipIncrements()
{
  Real temp = _temp[_qp];
  
  // Critical resolved shear stress in the input file
  // refers always to room temperature
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _gssT[i] = ( _dCRSS_dT_A + _dCRSS_dT_B * std::exp(- _dCRSS_dT_C * (temp - 293.0))) * 
	           _gss_tmp[i];
  }
  
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _slip_incr(i) = _a0(i) * 
	                std::pow(std::abs(_tau(i) / _gssT[i]), 1.0 / _xm(i)) *
                    std::copysign(1.0, _tau(i)) * _dt;
    if (std::abs(_slip_incr(i)) > _slip_incr_tol)
    {
      _err_tol = true;
#ifdef DEBUG
      mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
#endif
      return;
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
    _dslipdtau(i) = _a0(i) / _xm(i) *
                    std::pow(std::abs(_tau(i) / _gssT[i]), 1.0 / _xm(i) - 1.0) / _gssT[i] *
                    _dt;
}

