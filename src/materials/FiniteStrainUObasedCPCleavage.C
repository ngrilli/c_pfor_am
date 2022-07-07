// Nicolò Grilli
// University of Bristol
// 15 Febbraio 2022

#include "FiniteStrainUObasedCPCleavage.h"

#include "petscblaslapack.h"
#include "MooseException.h"
#include "MathUtils.h"

registerMooseObject("TensorMechanicsApp", FiniteStrainUObasedCPCleavage);

InputParameters
FiniteStrainUObasedCPCleavage::validParams()
{
  InputParameters params = FiniteStrainUObasedCP::validParams();
  params.addClassDescription("UserObject based Crystal Plasticity system with damage. "
                             "Free energy is decomposed into volumetric and non-volumetric " 
							 "parts as in: "
							 "Nicolo Grilli and Marisol Koslowski "
							 "The effect of crystal anisotropy and plastic response "
							 "on the dynamic fracture of energetic materials "
							 "Journal of Applied Physics 126, 155101 (2019).");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>(
      "E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<Real>("bulk_modulus_ref",0.0,"reference bulk modulus for vol/non-vol decomposition");
  return params;
}

FiniteStrainUObasedCPCleavage::FiniteStrainUObasedCPCleavage(const InputParameters & parameters)
  : FiniteStrainUObasedCP(parameters),
    _c(coupledValue("c")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _H(declareProperty<Real>("hist")), // History variable to avoid damage decrease 
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _E(declareProperty<Real>(getParam<MaterialPropertyName>("E_name"))),
    _dEdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"),
                                          getVar("c", 0)->name())),
    _d2Ed2c(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _dstress_dc(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())), 
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),		
	_D(getMaterialProperty<Real>("D_name")),
    _dDdc(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name())),
    _d2Dd2c(getMaterialPropertyDerivative<Real>(
        "D_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _bulk_modulus_ref(getParam<Real>("bulk_modulus_ref")), // reference bulk modulus for vol/non-vol decomposition
    _slip_plane_normals(declareProperty<std::vector<Real>>("slip_plane_normals")) // Slip plane normals
{
  _err_tol = false;

  _delta_dfgrd.zero();

  // resize the material properties for each userobject
  _mat_prop_slip_rates.resize(_num_uo_slip_rates);
  _mat_prop_slip_resistances.resize(_num_uo_slip_resistances);
  _mat_prop_state_vars.resize(_num_uo_state_vars);
  _mat_prop_state_vars_old.resize(_num_uo_state_vars);
  _mat_prop_state_var_evol_rate_comps.resize(_num_uo_state_var_evol_rate_comps);

  // resize the flow direction
  _flow_direction.resize(_num_uo_slip_rates);

  // resize local state variables
  _state_vars_old.resize(_num_uo_state_vars);
  _state_vars_old_stored.resize(_num_uo_state_vars);
  _state_vars_prev.resize(_num_uo_state_vars);

  // resize user objects
  _uo_slip_rates.resize(_num_uo_slip_rates);
  _uo_slip_resistances.resize(_num_uo_slip_resistances);
  _uo_state_vars.resize(_num_uo_state_vars);
  _uo_state_var_evol_rate_comps.resize(_num_uo_state_var_evol_rate_comps);

  // assign the user objects
  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
  {
    _uo_slip_rates[i] = &getUserObjectByName<CrystalPlasticitySlipRateCleavage>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_rates")[i]);
    _mat_prop_slip_rates[i] = &declareProperty<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_rates")[i]);
    _flow_direction[i] = &declareProperty<std::vector<RankTwoTensor>>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_rates")[i] + "_flow_direction");
  }

  for (unsigned int i = 0; i < _num_uo_slip_resistances; ++i)
  {
    _uo_slip_resistances[i] = &getUserObjectByName<CrystalPlasticitySlipResistance>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_resistances")[i]);
    _mat_prop_slip_resistances[i] = &declareProperty<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_resistances")[i]);
  }

  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
  {
    _uo_state_vars[i] = &getUserObjectByName<CrystalPlasticityStateVariable>(
        parameters.get<std::vector<UserObjectName>>("uo_state_vars")[i]);
    _mat_prop_state_vars[i] = &declareProperty<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_state_vars")[i]);
    _mat_prop_state_vars_old[i] = &getMaterialPropertyOld<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_state_vars")[i]);
  }

  for (unsigned int i = 0; i < _num_uo_state_var_evol_rate_comps; ++i)
  {
    _uo_state_var_evol_rate_comps[i] = &getUserObjectByName<CrystalPlasticityStateVarRateComponent>(
        parameters.get<std::vector<UserObjectName>>("uo_state_var_evol_rate_comps")[i]);
    _mat_prop_state_var_evol_rate_comps[i] = &declareProperty<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_state_var_evol_rate_comps")[i]);
  }

  _substep_dt = 0.0;
}

/**
 * Solves stress residual equation using NR.
 * Updates slip system resistances iteratively.
 * calcFlowDirection is modified to output _slip_plane_normals
 * to ACInterfaceSlipPlaneFracture for cleavage
 * along the slip plane
 */
void
FiniteStrainUObasedCPCleavage::computeQpStress()
{
  // number of slip systems
  Real NSlipSys;
  
  // Temporary vector to store slip plane normals
  // before assigning to the material property
  std::vector<Real> slip_plane_normals;
	
  // Userobject based crystal plasticity does not support face/boundary material property
  // calculation.
  if (isBoundaryMaterial())
    return;
  // Depth of substepping; Limited to maximum substep iteration
  unsigned int substep_iter = 1;
  // Calculated from substep_iter as 2^substep_iter
  unsigned int num_substep = 1;

  _dfgrd_tmp_old = _deformation_gradient_old[_qp];
  if (_dfgrd_tmp_old.det() == 0)
    _dfgrd_tmp_old.addIa(1.0);

  _delta_dfgrd = _deformation_gradient[_qp] - _dfgrd_tmp_old;

  // Saves the old stateful properties that are modified during sub stepping
  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
    _state_vars_old[i] = (*_mat_prop_state_vars_old[i])[_qp];

  // Set the size of _slip_plane_normals
  // variableSize() gives the number of slip systems
  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i) {
    NSlipSys = _uo_slip_rates[i]->variableSize();
  }
  _slip_plane_normals[_qp].resize(NSlipSys * LIBMESH_DIM);
  slip_plane_normals.resize(NSlipSys * LIBMESH_DIM);
  
  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
    _uo_slip_rates[i]->calcFlowDirection(_qp, (*_flow_direction[i])[_qp], slip_plane_normals);

  _slip_plane_normals[_qp] = slip_plane_normals;

  do
  {
    _err_tol = false;

    preSolveQp();

    _substep_dt = _dt / num_substep;

    for (unsigned int istep = 0; istep < num_substep; ++istep)
    {
      _dfgrd_tmp = (static_cast<Real>(istep) + 1) / num_substep * _delta_dfgrd + _dfgrd_tmp_old;

      solveQp();

      if (_err_tol)
      {
        substep_iter++;
        num_substep *= 2;
        break;
      }
    }

    if (substep_iter > _max_substep_iter && _err_tol)
      throw MooseException("FiniteStrainUObasedCP: Constitutive failure.");
  } while (_err_tol);

  postSolveQp();
}

void
FiniteStrainUObasedCPCleavage::calcResidual()
{
  RankTwoTensor iden(RankTwoTensor::initIdentity), ce, ee, ce_pk2, eqv_slip_incr, pk2_new;
  
  Real F_pos, F_neg; // tensile and compressive part of the elastic strain energy

  getSlipRates();
  if (_err_tol)
    return;

  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
    for (unsigned int j = 0; j < _uo_slip_rates[i]->variableSize(); ++j)
      eqv_slip_incr +=
          (*_flow_direction[i])[_qp][j] * (*_mat_prop_slip_rates[i])[_qp][j] * _substep_dt;

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;

  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;
  
  // Decompose ee into volumetric and non-volumetric
  // and calculate elastic energy and stress
  computeStrainVolumetric(F_pos, F_neg, ee, ce, pk2_new);
  
  // calculate history variable and
  // assign elastic free energy to _E
  // for the fracture model
  computeHistoryVariable(F_pos, F_neg);

  // Anisotropic undamaged
  // pk2_new = _elasticity_tensor[_qp] * ee;

  _resid = _pk2[_qp] - pk2_new;
}

void
FiniteStrainUObasedCPCleavage::computeStrainVolumetric(Real & F_pos, Real & F_neg, 
                                                        RankTwoTensor & ee, RankTwoTensor & ce, 
														RankTwoTensor & pk2_new)
{
  //Anisotropic elasticity
  Real Kb = 0.0; // Kb is the reference bulk modulus
  Real Je; // Je is relative elastic volume change
  Real Je23; // Je^{2/3}
  Real delta; // delta is the trace of volumetric part of elastic Green-Lagrange strain
  
  // Positive and negative volumetric free energies
  Real a_pos_vol, a_neg_vol;
  
  // Isochoric free energy
  Real a_pos_cpl;
  
  // Undamaged elastic energy tensor
  RankTwoTensor undamaged_elastic_energy;
  
  // Positive and negative part of the Piola-Kirchhoff stress
  RankTwoTensor pk2_pos;
  RankTwoTensor pk2_neg;
  
  // inverse of the right Cauchy–Green deformation tensor (elastic)
  RankTwoTensor invce;
  
  // Anisotropic elasticity (Luscher 2017)
  // Kb = K in (Luscher 2017)
  // Kb = (1/9) I : C : I
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      Kb += _elasticity_tensor[_qp](i, i, j, j);
  
  Kb = Kb / 9.0;

  Je = _fe.det();
  Je23 = std::pow(Je,2.0/3.0);
  delta = 1.5 * (Je23 - 1.0);
  
  // Calculate volumetric free energy
  // Equations 15 and 16 in Grilli, Koslowski, 2019
  
  if (Je >= 1.0) { // expansion
	  
    a_pos_vol = 0.5 * Kb * delta * delta;
    a_neg_vol = 0.0;
	  
  } else { // compression
	  
    a_pos_vol = 0.0;
	a_neg_vol = 0.5 * Kb * delta * delta;
	
  }
  
  // Calculate isochoric free energy
  // Equation 17 in Grilli, Koslowski, 2019

  undamaged_elastic_energy = _elasticity_tensor[_qp] * ee;
  undamaged_elastic_energy = 0.5 * ee * undamaged_elastic_energy;

  a_pos_cpl = undamaged_elastic_energy.trace();
  a_pos_cpl -= 0.5 * Kb * delta * delta;

  // Calculate positive and negative parts of the Piola-Kirchhoff stress
  // Equation 18 in Grilli, Koslowski, 2019
  if (Je >= 1.0) { // expansion
  
    pk2_pos = _elasticity_tensor[_qp] * ee;
    pk2_neg = 0.0;
	
  } else { // compression
  
    invce = ce.inverse();
	
	pk2_neg = Je23 * Kb * delta * invce;
	
    pk2_pos = (-1.0) * pk2_neg;	
    pk2_pos += _elasticity_tensor[_qp] * ee;

  }

  // Positive part of the stress is degraded by _D[_qp]
  pk2_new = (_D[_qp] * pk2_pos) + pk2_neg;
  
  // Positive and negative parts of the free energy
  // Equations 13 and 14 in Grilli, Koslowski, 2019
  F_pos = a_pos_vol + a_pos_cpl;
  F_neg = a_neg_vol;
  
  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = pk2_pos * _dDdc[_qp];
  
  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history variable
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = pk2_pos * _dDdc[_qp];
  
  // _Jacobian_mult is already defined in the CP base class
  
}

// compute history variable and assign to _E
// which is used by the fracture model for damage growth
// Damage grows only because of the positive part of the elastic energy F_pos
void
FiniteStrainUObasedCPCleavage::computeHistoryVariable(Real & F_pos, Real & F_neg)
{
  // Assign history variable
  Real hist_variable = _H_old[_qp];
  
  // _use_snes_vi_solver option not implemented

  if (F_pos > _H_old[_qp])
    _H[_qp] = F_pos;
  else
    _H[_qp] = _H_old[_qp];

  if (_use_current_hist)
    hist_variable = _H[_qp];

  // _barrier not implemented

  // Elastic free energy density and derivatives
  _E[_qp] = hist_variable * _D[_qp] + F_neg;
  _dEdc[_qp] = hist_variable * _dDdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp];

}

// update jacobian_mult by taking into account of the exact elasto-plastic tangent moduli
// it includes damage
void
FiniteStrainUObasedCPCleavage::elastoPlasticTangentModuli()
{
  RankFourTensor tan_mod;
  RankTwoTensor pk2fet, fepk2;
  RankFourTensor deedfe, dsigdpk2dfe, dfedf;

  // Fill in the matrix stiffness material property
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _fe(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _fe(k, i) * 0.5;
      }
	  
  usingTensorIndices(i_, j_, k_, l_);
  // This equation is exact for Je >= 1 but only an approximation for Je < 1
  dsigdpk2dfe = _fe.times<i_, k_, j_, l_>(_fe) * _D[_qp] * _elasticity_tensor[_qp] * deedfe;

  pk2fet = _pk2[_qp] * _fe.transpose();
  fepk2 = _fe * _pk2[_qp];

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
      {
        tan_mod(i, j, i, l) += pk2fet(l, j);
        tan_mod(i, j, j, l) += fepk2(i, l);
      }

  tan_mod += dsigdpk2dfe;

  Real je = _fe.det();
  if (je > 0.0)
    tan_mod /= je;

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
        dfedf(i, j, i, l) = _fp_inv(l, j);

  _Jacobian_mult[_qp] = tan_mod * dfedf;
}

// update jacobian_mult
// These are approximated tangent moduli
// but damage is included to make it more precise
void
FiniteStrainUObasedCPCleavage::elasticTangentModuli()
{
  _Jacobian_mult[_qp] = _D[_qp] * _elasticity_tensor[_qp];
}


