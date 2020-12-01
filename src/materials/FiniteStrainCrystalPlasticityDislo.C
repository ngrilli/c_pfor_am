// Nicolo Grilli
// Daijun Hu 
// National University of Singapore
// 16 Novembre 2020

#include "FiniteStrainCrystalPlasticityDislo.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

#include <fstream>
#include <cmath>

registerMooseObject("TensorMechanicsApp", FiniteStrainCrystalPlasticityDislo);

InputParameters
FiniteStrainCrystalPlasticityDislo::validParams()
{
  InputParameters params = FiniteStrainCrystalPlasticity::validParams();
  params.addClassDescription("Crystal Plasticity with thermal eigenstrain. "
                             "Temperature dependence of the CRSS. "
							 "Dislocation based model. "
							 "Stress dependent dislocation velocity. ");
  params.addCoupledVar("temp",293.0,"Temperature");
  params.addCoupledVar("rho_edge_pos_1",0.0,"Positive edge dislocation density: slip system 1");
  params.addCoupledVar("rho_edge_neg_1",0.0,"Negative edge dislocation density: slip system 1");
  params.addCoupledVar("rho_edge_pos_2",0.0,"Positive edge dislocation density: slip system 2");
  params.addCoupledVar("rho_edge_neg_2",0.0,"Negative edge dislocation density: slip system 2");
  params.addCoupledVar("rho_edge_pos_3",0.0,"Positive edge dislocation density: slip system 3");
  params.addCoupledVar("rho_edge_neg_3",0.0,"Negative edge dislocation density: slip system 3");
  params.addCoupledVar("rho_edge_pos_4",0.0,"Positive edge dislocation density: slip system 4");
  params.addCoupledVar("rho_edge_neg_4",0.0,"Negative edge dislocation density: slip system 4");
  params.addCoupledVar("rho_edge_pos_5",0.0,"Positive edge dislocation density: slip system 5");
  params.addCoupledVar("rho_edge_neg_5",0.0,"Negative edge dislocation density: slip system 5");
  params.addCoupledVar("rho_edge_pos_6",0.0,"Positive edge dislocation density: slip system 6");
  params.addCoupledVar("rho_edge_neg_6",0.0,"Negative edge dislocation density: slip system 6");
  params.addCoupledVar("rho_edge_pos_7",0.0,"Positive edge dislocation density: slip system 7");
  params.addCoupledVar("rho_edge_neg_7",0.0,"Negative edge dislocation density: slip system 7");
  params.addCoupledVar("rho_edge_pos_8",0.0,"Positive edge dislocation density: slip system 8");
  params.addCoupledVar("rho_edge_neg_8",0.0,"Negative edge dislocation density: slip system 8");
  params.addCoupledVar("rho_edge_pos_9",0.0,"Positive edge dislocation density: slip system 9");
  params.addCoupledVar("rho_edge_neg_9",0.0,"Negative edge dislocation density: slip system 9");
  params.addCoupledVar("rho_edge_pos_10",0.0,"Positive edge dislocation density: slip system 10");
  params.addCoupledVar("rho_edge_neg_10",0.0,"Negative edge dislocation density: slip system 10");
  params.addCoupledVar("rho_edge_pos_11",0.0,"Positive edge dislocation density: slip system 11");
  params.addCoupledVar("rho_edge_neg_11",0.0,"Negative edge dislocation density: slip system 11");
  params.addCoupledVar("rho_edge_pos_12",0.0,"Positive edge dislocation density: slip system 12");
  params.addCoupledVar("rho_edge_neg_12",0.0,"Negative edge dislocation density: slip system 12");
  params.addCoupledVar("rho_screw_pos_1",0.0,"Positive screw dislocation density: slip system 1");
  params.addCoupledVar("rho_screw_neg_1",0.0,"Negative screw dislocation density: slip system 1");
  params.addCoupledVar("rho_screw_pos_2",0.0,"Positive screw dislocation density: slip system 2");
  params.addCoupledVar("rho_screw_neg_2",0.0,"Negative screw dislocation density: slip system 2");
  params.addCoupledVar("rho_screw_pos_3",0.0,"Positive screw dislocation density: slip system 3");
  params.addCoupledVar("rho_screw_neg_3",0.0,"Negative screw dislocation density: slip system 3");
  params.addCoupledVar("rho_screw_pos_4",0.0,"Positive screw dislocation density: slip system 4");
  params.addCoupledVar("rho_screw_neg_4",0.0,"Negative screw dislocation density: slip system 4");
  params.addCoupledVar("rho_screw_pos_5",0.0,"Positive screw dislocation density: slip system 5");
  params.addCoupledVar("rho_screw_neg_5",0.0,"Negative screw dislocation density: slip system 5");
  params.addCoupledVar("rho_screw_pos_6",0.0,"Positive screw dislocation density: slip system 6");
  params.addCoupledVar("rho_screw_neg_6",0.0,"Negative screw dislocation density: slip system 6");
  params.addCoupledVar("rho_screw_pos_7",0.0,"Positive screw dislocation density: slip system 7");
  params.addCoupledVar("rho_screw_neg_7",0.0,"Negative screw dislocation density: slip system 7");
  params.addCoupledVar("rho_screw_pos_8",0.0,"Positive screw dislocation density: slip system 8");
  params.addCoupledVar("rho_screw_neg_8",0.0,"Negative screw dislocation density: slip system 8");
  params.addCoupledVar("rho_screw_pos_9",0.0,"Positive screw dislocation density: slip system 9");
  params.addCoupledVar("rho_screw_neg_9",0.0,"Negative screw dislocation density: slip system 9");
  params.addCoupledVar("rho_screw_pos_10",0.0,"Positive screw dislocation density: slip system 10");
  params.addCoupledVar("rho_screw_neg_10",0.0,"Negative screw dislocation density: slip system 10");
  params.addCoupledVar("rho_screw_pos_11",0.0,"Positive screw dislocation density: slip system 11");
  params.addCoupledVar("rho_screw_neg_11",0.0,"Negative screw dislocation density: slip system 11");
  params.addCoupledVar("rho_screw_pos_12",0.0,"Positive screw dislocation density: slip system 12");
  params.addCoupledVar("rho_screw_neg_12",0.0,"Negative screw dislocation density: slip system 12");
  params.addParam<Real>("thermal_expansion",0.0,"Thermal expansion coefficient");
  params.addParam<Real>("reference_temperature",293.0,"reference temperature for thermal expansion");
  params.addParam<Real>("dCRSS_dT_A",1.0,"A coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 293.0))");
  params.addParam<Real>("dCRSS_dT_B",0.0,"B coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 293.0))");
  params.addParam<Real>("dCRSS_dT_C",0.0,"C coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 293.0))");
  params.addParam<Real>("dislo_mobility",0.0,"Dislocation mobility");
  params.addParam<Real>("burgers_vector_mag",0.0,"Magnitude of the Burgers vector");
  params.addParam<Real>("shear_modulus_hardening",86000.0,"Shear modulus in Taylor hardening law");
  return params;
}

FiniteStrainCrystalPlasticityDislo::FiniteStrainCrystalPlasticityDislo(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
    _temp(coupledValue("temp")),
    _rho_edge_pos_1(coupledValue("rho_edge_pos_1")),
	_rho_edge_neg_1(coupledValue("rho_edge_neg_1")),
    _rho_edge_pos_2(coupledValue("rho_edge_pos_2")),
	_rho_edge_neg_2(coupledValue("rho_edge_neg_2")),
    _rho_edge_pos_3(coupledValue("rho_edge_pos_3")),
	_rho_edge_neg_3(coupledValue("rho_edge_neg_3")),
    _rho_edge_pos_4(coupledValue("rho_edge_pos_4")),
	_rho_edge_neg_4(coupledValue("rho_edge_neg_4")),
    _rho_edge_pos_5(coupledValue("rho_edge_pos_5")),
	_rho_edge_neg_5(coupledValue("rho_edge_neg_5")),
    _rho_edge_pos_6(coupledValue("rho_edge_pos_6")),
	_rho_edge_neg_6(coupledValue("rho_edge_neg_6")),
    _rho_edge_pos_7(coupledValue("rho_edge_pos_7")),
	_rho_edge_neg_7(coupledValue("rho_edge_neg_7")),
    _rho_edge_pos_8(coupledValue("rho_edge_pos_8")),
	_rho_edge_neg_8(coupledValue("rho_edge_neg_8")),
    _rho_edge_pos_9(coupledValue("rho_edge_pos_9")),
	_rho_edge_neg_9(coupledValue("rho_edge_neg_9")),
    _rho_edge_pos_10(coupledValue("rho_edge_pos_10")),
	_rho_edge_neg_10(coupledValue("rho_edge_neg_10")),
    _rho_edge_pos_11(coupledValue("rho_edge_pos_11")),
	_rho_edge_neg_11(coupledValue("rho_edge_neg_11")),
    _rho_edge_pos_12(coupledValue("rho_edge_pos_12")),
    _rho_edge_neg_12(coupledValue("rho_edge_neg_12")),	
    _rho_screw_pos_1(coupledValue("rho_screw_pos_1")),
	_rho_screw_neg_1(coupledValue("rho_screw_neg_1")),
    _rho_screw_pos_2(coupledValue("rho_screw_pos_2")),
	_rho_screw_neg_2(coupledValue("rho_screw_neg_2")),
    _rho_screw_pos_3(coupledValue("rho_screw_pos_3")),
	_rho_screw_neg_3(coupledValue("rho_screw_neg_3")),
    _rho_screw_pos_4(coupledValue("rho_screw_pos_4")),
	_rho_screw_neg_4(coupledValue("rho_screw_neg_4")),
    _rho_screw_pos_5(coupledValue("rho_screw_pos_5")),
	_rho_screw_neg_5(coupledValue("rho_screw_neg_5")),
    _rho_screw_pos_6(coupledValue("rho_screw_pos_6")),
	_rho_screw_neg_6(coupledValue("rho_screw_neg_6")),
    _rho_screw_pos_7(coupledValue("rho_screw_pos_7")),
	_rho_screw_neg_7(coupledValue("rho_screw_neg_7")),
    _rho_screw_pos_8(coupledValue("rho_screw_pos_8")),
	_rho_screw_neg_8(coupledValue("rho_screw_neg_8")),
    _rho_screw_pos_9(coupledValue("rho_screw_pos_9")),
	_rho_screw_neg_9(coupledValue("rho_screw_neg_9")),
    _rho_screw_pos_10(coupledValue("rho_screw_pos_10")),
	_rho_screw_neg_10(coupledValue("rho_screw_neg_10")),
    _rho_screw_pos_11(coupledValue("rho_screw_pos_11")),
	_rho_screw_neg_11(coupledValue("rho_screw_neg_11")),
    _rho_screw_pos_12(coupledValue("rho_screw_pos_12")),
    _rho_screw_neg_12(coupledValue("rho_screw_neg_12")),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _dCRSS_dT_A(getParam<Real>("dCRSS_dT_A")),
	_dCRSS_dT_B(getParam<Real>("dCRSS_dT_B")),
	_dCRSS_dT_C(getParam<Real>("dCRSS_dT_C")),
	_dislo_mobility(getParam<Real>("dislo_mobility")),
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")), // Magnitude of the Burgers vector
	_shear_modulus_hardening(getParam<Real>("shear_modulus_hardening")), // Shear modulus in Taylor hardening law
	_gssT(_nss),
    _edge_slip_direction(declareProperty<std::vector<Real>>("edge_slip_direction")), // Edge slip directions
	_screw_slip_direction(declareProperty<std::vector<Real>>("screw_slip_direction")), // Screw slip direction
	_slip_incr_out(declareProperty<std::vector<Real>>("slip_incr_out")), // Slip system resistances
	_dislo_velocity(declareProperty<std::vector<Real>>("dislo_velocity")), // Dislocation velocity
	_ddislo_velocity_dtau(declareProperty<std::vector<Real>>("ddislo_velocity_dtau")) // Derivative of dislo velocity
{	
}

void
FiniteStrainCrystalPlasticityDislo::calcResidual( RankTwoTensor &resid )
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

  // Introduce temperature dependence of the CRSS
  TempDependCRSS();

  // calculate dislocation velocity
  // and store it for advection kernel
  // necessary to call it here because getSlipIncrements
  // depends on dislocation velocity, changing at each iteration
  // of the CP algorithm
  getDisloVelocity();

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
  
  // It would be better to call this function in postSolveQp()
  // so it is not called more times than necessary
  OutputSlipDirection();
}

// Critical resolved shear stress decreases exponentially with temperature
// A + B exp(- C * (T - 293.0))
void
FiniteStrainCrystalPlasticityDislo::TempDependCRSS()
{
  Real temp = _temp[_qp];
  
  // Critical resolved shear stress in the input file
  // refers always to room temperature
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _gssT[i] = ( _dCRSS_dT_A + _dCRSS_dT_B * std::exp(- _dCRSS_dT_C * (temp - 293.0))) * 
	           _gss_tmp[i];
  }
}

// Calculate slip increment,dslipdtau
void
FiniteStrainCrystalPlasticityDislo::getSlipIncrements()
{  
  std::vector<Real> rho_edge_pos(_nss);
  std::vector<Real> rho_edge_neg(_nss);
  std::vector<Real> rho_screw_pos(_nss);
  std::vector<Real> rho_screw_neg(_nss);

  // Assign dislocation density vectors
  rho_edge_pos[0] = _rho_edge_pos_1[_qp];
  rho_edge_pos[1] = _rho_edge_pos_2[_qp];
  rho_edge_pos[2] = _rho_edge_pos_3[_qp];
  rho_edge_pos[3] = _rho_edge_pos_4[_qp];
  rho_edge_pos[4] = _rho_edge_pos_5[_qp];
  rho_edge_pos[5] = _rho_edge_pos_6[_qp];
  rho_edge_pos[6] = _rho_edge_pos_7[_qp];
  rho_edge_pos[7] = _rho_edge_pos_8[_qp];
  rho_edge_pos[8] = _rho_edge_pos_9[_qp];
  rho_edge_pos[9] = _rho_edge_pos_10[_qp];
  rho_edge_pos[10] = _rho_edge_pos_11[_qp];
  rho_edge_pos[11] = _rho_edge_pos_12[_qp];
  
  rho_edge_neg[0] = _rho_edge_neg_1[_qp];
  rho_edge_neg[1] = _rho_edge_neg_2[_qp];
  rho_edge_neg[2] = _rho_edge_neg_3[_qp];
  rho_edge_neg[3] = _rho_edge_neg_4[_qp];
  rho_edge_neg[4] = _rho_edge_neg_5[_qp];
  rho_edge_neg[5] = _rho_edge_neg_6[_qp];
  rho_edge_neg[6] = _rho_edge_neg_7[_qp];
  rho_edge_neg[7] = _rho_edge_neg_8[_qp];
  rho_edge_neg[8] = _rho_edge_neg_9[_qp];
  rho_edge_neg[9] = _rho_edge_neg_10[_qp];
  rho_edge_neg[10] = _rho_edge_neg_11[_qp];
  rho_edge_neg[11] = _rho_edge_neg_12[_qp];
  
  rho_screw_pos[0] = _rho_screw_pos_1[_qp];
  rho_screw_pos[1] = _rho_screw_pos_2[_qp];
  rho_screw_pos[2] = _rho_screw_pos_3[_qp];
  rho_screw_pos[3] = _rho_screw_pos_4[_qp];
  rho_screw_pos[4] = _rho_screw_pos_5[_qp];
  rho_screw_pos[5] = _rho_screw_pos_6[_qp];
  rho_screw_pos[6] = _rho_screw_pos_7[_qp];
  rho_screw_pos[7] = _rho_screw_pos_8[_qp];
  rho_screw_pos[8] = _rho_screw_pos_9[_qp];
  rho_screw_pos[9] = _rho_screw_pos_10[_qp];
  rho_screw_pos[10] = _rho_screw_pos_11[_qp];
  rho_screw_pos[11] = _rho_screw_pos_12[_qp];
  
  rho_screw_neg[0] = _rho_screw_neg_1[_qp];
  rho_screw_neg[1] = _rho_screw_neg_2[_qp];
  rho_screw_neg[2] = _rho_screw_neg_3[_qp];
  rho_screw_neg[3] = _rho_screw_neg_4[_qp];
  rho_screw_neg[4] = _rho_screw_neg_5[_qp];
  rho_screw_neg[5] = _rho_screw_neg_6[_qp];
  rho_screw_neg[6] = _rho_screw_neg_7[_qp];
  rho_screw_neg[7] = _rho_screw_neg_8[_qp];
  rho_screw_neg[8] = _rho_screw_neg_9[_qp];
  rho_screw_neg[9] = _rho_screw_neg_10[_qp];
  rho_screw_neg[10] = _rho_screw_neg_11[_qp];
  rho_screw_neg[11] = _rho_screw_neg_12[_qp];

  // Positive and negative dislocation give the same
  // contribution to Lp even if their velocity is opposite
  // same for edge and screw
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _slip_incr(i) = (rho_edge_pos[i] + rho_edge_neg[i] + rho_screw_pos[i] + rho_screw_neg[i]) * 
	                std::abs(_dislo_velocity[_qp][i]) * _burgers_vector_mag *
                    std::copysign(1.0, _tau(i)) * _dt;
  }

  // Derivative is always positive
  for (unsigned int i = 0; i < _nss; ++i)
    _dslipdtau(i) = (rho_edge_pos[i] + rho_edge_neg[i] + rho_screw_pos[i] + rho_screw_neg[i]) * 
                    _ddislo_velocity_dtau[_qp][i] * _burgers_vector_mag * _dt;
		
  for (unsigned int i = 0; i < _nss; ++i)
  {
    if (std::abs(_slip_incr(i)) > _slip_incr_tol)
    {
      //_err_tol = true;
      //mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
	  
	  _slip_incr(i) = _slip_incr_tol * std::copysign(1.0, _tau(i));
	  
    }	  
  }	  
					
  // store slip increment for output
  _slip_incr_out[_qp].resize(_nss);
  
  for (unsigned int i = 0; i < _nss; ++i) {
    _slip_incr_out[_qp][i] = _slip_incr(i);
  }
}

// Calculate dislocation velocity (edge and screw) as a function
// of the resolved shear stress and its derivative
void
FiniteStrainCrystalPlasticityDislo::getDisloVelocity()
{
  _dislo_velocity[_qp].resize(_nss);
  _ddislo_velocity_dtau[_qp].resize(_nss);
  
  for (unsigned int i = 0; i < _nss; ++i)
  {
	_dislo_velocity[_qp][i] = 0.0;
	_ddislo_velocity_dtau[_qp][i] = 0.0;
	
	if (std::abs(_tau(i)) > _gssT[i]) {
      _dislo_velocity[_qp][i] = _dislo_mobility * (std::abs(_tau(i)) - _gssT[i])
	                            * std::copysign(1.0, _tau(i));
	  // Derivative is always positive
	  _ddislo_velocity_dtau[_qp][i] = _dislo_mobility;
	}
  }
}

// Store slip direction
// to couple with dislocation transport
void
FiniteStrainCrystalPlasticityDislo::OutputSlipDirection()
{
  DenseVector<Real> mo(LIBMESH_DIM * _nss);
  DenseVector<Real> no(LIBMESH_DIM * _nss);
  
  // Temporary directions and normals to calculate
  // screw dislocation slip direction
  RealVectorValue temp_mo;
  RealVectorValue temp_no;
  RealVectorValue temp_screw_mo;

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
	
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        no(i * LIBMESH_DIM + j) =
            no(i * LIBMESH_DIM + j) + _crysrot[_qp](j, k) * _no(i * LIBMESH_DIM + k);
    }
  }
 
  _edge_slip_direction[_qp].resize(LIBMESH_DIM * _nss);
  _screw_slip_direction[_qp].resize(LIBMESH_DIM * _nss);

  // Store slip direction (already normalized)
  // for edge and screw dislocations
  // to couple with dislocation transport
  for (unsigned int i = 0; i < _nss; ++i)
  {
	for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	{
	  temp_mo(j) = mo(i * LIBMESH_DIM + j);
	  temp_no(j) = no(i * LIBMESH_DIM + j);
	}		
	
	temp_screw_mo = temp_mo.cross(temp_no);
	  
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
  	  _edge_slip_direction[_qp][i * LIBMESH_DIM + j] = mo(i * LIBMESH_DIM + j);
	  _screw_slip_direction[_qp][i * LIBMESH_DIM + j] = temp_screw_mo(j);
  	}
  }
  
}

/**
 * Calculate slip system resistance (CRSS)
 * based on Taylor hardening model
 */
void
FiniteStrainCrystalPlasticityDislo::updateGss()
{
  Real qab; // Taylor hardening
  Real TotalRho = 0.0; // total dislocation density
  
  std::vector<Real> rho_edge_pos(_nss);
  std::vector<Real> rho_edge_neg(_nss);
  std::vector<Real> rho_screw_pos(_nss);
  std::vector<Real> rho_screw_neg(_nss);

  // Assign dislocation density vectors
  rho_edge_pos[0] = _rho_edge_pos_1[_qp];
  rho_edge_pos[1] = _rho_edge_pos_2[_qp];
  rho_edge_pos[2] = _rho_edge_pos_3[_qp];
  rho_edge_pos[3] = _rho_edge_pos_4[_qp];
  rho_edge_pos[4] = _rho_edge_pos_5[_qp];
  rho_edge_pos[5] = _rho_edge_pos_6[_qp];
  rho_edge_pos[6] = _rho_edge_pos_7[_qp];
  rho_edge_pos[7] = _rho_edge_pos_8[_qp];
  rho_edge_pos[8] = _rho_edge_pos_9[_qp];
  rho_edge_pos[9] = _rho_edge_pos_10[_qp];
  rho_edge_pos[10] = _rho_edge_pos_11[_qp];
  rho_edge_pos[11] = _rho_edge_pos_12[_qp];
  
  rho_edge_neg[0] = _rho_edge_neg_1[_qp];
  rho_edge_neg[1] = _rho_edge_neg_2[_qp];
  rho_edge_neg[2] = _rho_edge_neg_3[_qp];
  rho_edge_neg[3] = _rho_edge_neg_4[_qp];
  rho_edge_neg[4] = _rho_edge_neg_5[_qp];
  rho_edge_neg[5] = _rho_edge_neg_6[_qp];
  rho_edge_neg[6] = _rho_edge_neg_7[_qp];
  rho_edge_neg[7] = _rho_edge_neg_8[_qp];
  rho_edge_neg[8] = _rho_edge_neg_9[_qp];
  rho_edge_neg[9] = _rho_edge_neg_10[_qp];
  rho_edge_neg[10] = _rho_edge_neg_11[_qp];
  rho_edge_neg[11] = _rho_edge_neg_12[_qp];
  
  rho_screw_pos[0] = _rho_screw_pos_1[_qp];
  rho_screw_pos[1] = _rho_screw_pos_2[_qp];
  rho_screw_pos[2] = _rho_screw_pos_3[_qp];
  rho_screw_pos[3] = _rho_screw_pos_4[_qp];
  rho_screw_pos[4] = _rho_screw_pos_5[_qp];
  rho_screw_pos[5] = _rho_screw_pos_6[_qp];
  rho_screw_pos[6] = _rho_screw_pos_7[_qp];
  rho_screw_pos[7] = _rho_screw_pos_8[_qp];
  rho_screw_pos[8] = _rho_screw_pos_9[_qp];
  rho_screw_pos[9] = _rho_screw_pos_10[_qp];
  rho_screw_pos[10] = _rho_screw_pos_11[_qp];
  rho_screw_pos[11] = _rho_screw_pos_12[_qp];
  
  rho_screw_neg[0] = _rho_screw_neg_1[_qp];
  rho_screw_neg[1] = _rho_screw_neg_2[_qp];
  rho_screw_neg[2] = _rho_screw_neg_3[_qp];
  rho_screw_neg[3] = _rho_screw_neg_4[_qp];
  rho_screw_neg[4] = _rho_screw_neg_5[_qp];
  rho_screw_neg[5] = _rho_screw_neg_6[_qp];
  rho_screw_neg[6] = _rho_screw_neg_7[_qp];
  rho_screw_neg[7] = _rho_screw_neg_8[_qp];
  rho_screw_neg[8] = _rho_screw_neg_9[_qp];
  rho_screw_neg[9] = _rho_screw_neg_10[_qp];
  rho_screw_neg[10] = _rho_screw_neg_11[_qp];
  rho_screw_neg[11] = _rho_screw_neg_12[_qp];

  // Is this update necessary
  // for the constitutive model
  _accslip_tmp = _accslip_tmp_old;
  for (unsigned int i = 0; i < _nss; ++i)
    _accslip_tmp += std::abs(_slip_incr(i));

  for (unsigned int i = 0; i < _nss; ++i)
    TotalRho += (rho_edge_pos[i] + rho_edge_neg[i] + rho_screw_pos[i] + rho_screw_neg[i]);

  qab = 0.4 * _shear_modulus_hardening * _burgers_vector_mag * std::sqrt(TotalRho);
  
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _gss_tmp[i] = qab;
  }
}

