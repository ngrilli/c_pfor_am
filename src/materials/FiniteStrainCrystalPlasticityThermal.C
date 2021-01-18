// Nicolo Grilli
// Daijun Hu 
// National University of Singapore
// 27 Ottobre 2020

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
  params.addCoupledVar("temp",303.0,"Temperature");
  params.addParam<Real>("thermal_expansion",0.0,"Thermal expansion coefficient");
  params.addParam<Real>("reference_temperature",303.0,"reference temperature for thermal expansion");
  params.addParam<Real>("dCRSS_dT_A",1.0,"A coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 303.0))");
  params.addParam<Real>("dCRSS_dT_B",0.0,"B coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 303.0))");
  params.addParam<Real>("dCRSS_dT_C",0.0,"C coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 303.0))");
  params.addParam<Real>("dCTE_dT",0.0,"coefficient for the increase of thermal expansion coefficient");
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
    _dCTE_dT(getParam<Real>("dCTE_dT")),
	_gssT(_nss),
    _lattice_strain(declareProperty<RankTwoTensor>("lattice_strain")),
    _slip_direction(declareProperty<std::vector<Real>>("slip_direction")), // Slip directions
	_slip_incr_out(declareProperty<std::vector<Real>>("slip_incr_out"))    // Slip system resistances
{	
}

void
FiniteStrainCrystalPlasticityThermal::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;
  Real temp = _temp[_qp];
  Real thermal_expansion = _thermal_expansion; 
  Real reference_temperature = _reference_temperature;
  Real dCTE_dT =_dCTE_dT;

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
  
// CTE:alpha(T)=dCTE_dT*T+b
// Unit of temperature is K, b is the CTE value at 0K estimated by the linear fitting parts from reference
  RankTwoTensor thermal_eigenstrain;
  thermal_eigenstrain = (1.0 / 2.0) * (std::exp((2.0/3.0) 
                      * ((1.0/2.0) * dCTE_dT * (temp-reference_temperature)*(temp+reference_temperature)
                      + thermal_expansion * (temp - reference_temperature))) - 1.0) * iden;
  pk2_new = _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);
  
  resid = _pk2_tmp - pk2_new;
  
  // It would be better to call the following lines in postSolveQp()
  // so it is not called more times than necessary
  OutputSlipDirection();
  _lattice_strain[_qp] = _crysrot[_qp].transpose() * ee * _crysrot[_qp];
}

// Calculate slip increment,dslipdtau
// Critical resolved shear stress decreases exponentially with temperature
// A + B exp(- C * (T - 303.0))
void
FiniteStrainCrystalPlasticityThermal::getSlipIncrements()
{
  Real temp = _temp[_qp];
  
  // Critical resolved shear stress in the input file
  // refers always to room temperature
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _gssT[i] = ( _dCRSS_dT_A + _dCRSS_dT_B * std::exp(- _dCRSS_dT_C * (temp - _reference_temperature))) * 
	           _gss_tmp[i];
  }
  
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _slip_incr(i) = _a0(i) * 
	                std::pow(std::abs(_tau(i) / _gssT[i]), 1.0 / _xm(i)) *
                    std::copysign(1.0, _tau(i)) * _dt;
    if (std::abs(_slip_incr(i)) > _slip_incr_tol)
    {
      //_err_tol = true;
      //mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
	  
	  _slip_incr(i) = _slip_incr_tol * std::copysign(1.0, _tau(i));

    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dslipdtau(i) = _a0(i) / _xm(i) *
                    std::pow(std::abs(_tau(i) / _gssT[i]), 1.0 / _xm(i) - 1.0) / _gssT[i] *
                    _dt;
	if (std::abs(_slip_incr(i)) > 0.99*_slip_incr_tol) 
	{
	  _dslipdtau(i) = _a0(i) / _xm(i) *
					  std::pow(std::abs(_slip_incr_tol / _a0(i) / _dt), 1.0 - _xm(i)) / 
					  _gssT[i] * _dt;	
	}
  }

  // store slip increment for output
  _slip_incr_out[_qp].resize(_nss);
  
  for (unsigned int i = 0; i < _nss; ++i) {
    _slip_incr_out[_qp][i] = _slip_incr(i);
  }
  
}

// Store slip direction
// to couple with dislocation transport
void
FiniteStrainCrystalPlasticityThermal::OutputSlipDirection()
{
  DenseVector<Real> mo(LIBMESH_DIM * _nss);

  // Update slip direction with crystal orientation
  for (unsigned int i = 0; i < _nss; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo(i * LIBMESH_DIM + j) =
            mo(i * LIBMESH_DIM + j) + _crysrot[_qp](j, k) * _mo(i * LIBMESH_DIM + k);
    }
  }
 
  _slip_direction[_qp].resize(LIBMESH_DIM * _nss);

  // Store slip direction (already normalized)
  // to couple with dislocation transport
  for (unsigned int i = 0; i < _nss; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
  	  _slip_direction[_qp][i * LIBMESH_DIM + j] = mo(i * LIBMESH_DIM + j);
  	}
  }
  
}

