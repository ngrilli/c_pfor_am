// Nicolo Grilli
// University of Bristol
// 18 Ottobre 2021

// Constitutive model from:
// Nicolo Grilli , Alan C.F. Cocks , Edmund Tarleton
// Crystal plasticity finite element modelling of coarse-grained alpha-uranium
// Computational Materials Science 171 (2020) 109276
// Plus Armstrong-Frederick backstress term:
// Armstrong, P.J., Frederick, C.O., 1966. 
// A Mathematical Representation of the Multiaxial Bauschinger Effect, 
// G.E.G.B. Report RD/B/N. Central Electricity Generating Board.

#include "FiniteStrainCrystalPlasticityBackstress.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

#include <fstream>
#include <cmath>

registerMooseObject("TensorMechanicsApp", FiniteStrainCrystalPlasticityBackstress);

InputParameters
FiniteStrainCrystalPlasticityBackstress::validParams()
{
  InputParameters params = FiniteStrainCrystalPlasticity::validParams();
  params.addClassDescription("Crystal Plasticity with thermal eigenstrain"
							 "Dislocation densities are internal variables");
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
  params.addParam<Real>("k_bol", 1.38e-23, "Boltzmann constant"); 
  params.addParam<Real>("da0", 0.0, "Constant dislocation annihilation length");
  params.addParam<Real>("log_strain_rate_ratio", 1.0, "Logarithm of the strain rate ratio");
  params.addParam<Real>("drag_stress", 900.0, "Drag stress");
  params.addParam<Real>("ka", 0.0, "Pre-factor for dislocation multiplication");
  params.addParam<Real>("burgers_vector", 0.0, "Burgers vector magnitude");
  params.addParam<Real>("projected_mu", 74000.0, "Projected shear modulus on the slip systems");
  params.addParam<Real>("tau0", 7.0, "Constant friction stress");
  params.addParam<Real>("init_rho_for",1.0,"initial value of forest dislocation density, same values for all slip systems");
  params.addParam<Real>("init_rho_sub",1.0,"initial value of substructure dislocation density");
  params.addParam<bool>("rho_sub_flag",false,"Flag to determine whether to include rho_sub in simulations");
  params.addParam<Real>("c_backstress",0.0,"c parameter in Armstrong-Frederick backstress evolution");
  params.addParam<Real>("d_backstress",0.0,"d parameter in Armstrong-Frederick backstress evolution");
  return params;
}

FiniteStrainCrystalPlasticityBackstress::FiniteStrainCrystalPlasticityBackstress(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
    _temp(coupledValue("temp")),     
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _dCRSS_dT_A(getParam<Real>("dCRSS_dT_A")),
	_dCRSS_dT_B(getParam<Real>("dCRSS_dT_B")),
	_dCRSS_dT_C(getParam<Real>("dCRSS_dT_C")),
    _dCTE_dT(getParam<Real>("dCTE_dT")),
	_k_bol(getParam<Real>("k_bol")), // Boltzmann constant
	_da0(getParam<Real>("da0")), // Constant dislocation annihilation length
	_log_strain_rate_ratio(getParam<Real>("log_strain_rate_ratio")),
	_drag_stress(getParam<Real>("drag_stress")), // Drag stress
	_ka(getParam<Real>("ka")), // Pre-factor for dislocation multiplication
	_burgers_vector(getParam<Real>("burgers_vector")), // Burgers vector magnitude
	_projected_mu(getParam<Real>("projected_mu")), // Projected shear modulus on the slip systems
	_tau0(getParam<Real>("tau0")), // Constant friction stress
    _init_rho_for(getParam<Real>("init_rho_for")), // Initial value of forest dislocation density
    _init_rho_sub(getParam<Real>("init_rho_sub")), // Initial value of substructure dislocation density
    _rho_sub_flag(getParam<bool>("rho_sub_flag")), // Flag to determine whether to include rho_sub in the model
	_c_backstress(getParam<Real>("c_backstress")), // c parameter in Armstrong-Frederick backstress evolution
	_d_backstress(getParam<Real>("d_backstress")), // d parameter in Armstrong-Frederick backstress evolution
	_gssT(_nss),
    _lattice_strain(declareProperty<RankTwoTensor>("lattice_strain")),
	_slip_incr_out(declareProperty<std::vector<Real>>("slip_incr_out")),   // Slip system resistances
	_rho_for(declareProperty<std::vector<Real>>("rho_for")), // Forest dislocation density
	_rho_for_old(getMaterialPropertyOld<std::vector<Real>>("rho_for")), 
	_rho_sub(declareProperty<Real>("rho_sub")), // Substructure dislocation density
	_rho_sub_old(getMaterialPropertyOld<Real>("rho_sub")),
	_tau_b(declareProperty<std::vector<Real>>("tau_b")), // Backstress for each slip system
	_tau_b_old(getMaterialPropertyOld<std::vector<Real>>("tau_b")),
	_rho_for_tmp(_nss), // give dimension to temporary vector
	_rho_for_tmp_old(_nss), // give dimension to temporary vector
	_tau_b_tmp(_nss), // give dimension to temporary vector
	_tau_b_tmp_old(_nss) // give dimension to temporary vector
{	
}

// Initialize the state variables of the dislocation model
void
FiniteStrainCrystalPlasticityBackstress::initAdditionalProps()
{
  // give dimension to material property
  // old material property will take the same dimension automatically  
  _rho_for[_qp].resize(_nss);	
  _tau_b[_qp].resize(_nss);
	
  for (unsigned int i = 0; i < _nss; ++i) // initialise forest dislocation densities
    _rho_for[_qp][i] = _init_rho_for;
  
  _rho_sub[_qp] = _init_rho_sub;	// initialise substructure dislocation density
  
  for (unsigned int i = 0; i < _nss; ++i) // initialise backstress
    _tau_b[_qp][i] = 0.0;
  
}

// Assign state variables from the previous time step
// to temporary variables
void
FiniteStrainCrystalPlasticityBackstress::preSolveStatevar()
{
  if (_max_substep_iter == 1) // No substepping
  {
    _gss_tmp = _gss_old[_qp];
    _accslip_tmp_old = _acc_slip_old[_qp];
	
	for (unsigned int i = 0; i < _nss; ++i)
	  _rho_for_tmp_old[i] = _rho_for_old[_qp][i];
  
    _rho_sub_tmp_old = _rho_sub_old[_qp];
	
	for (unsigned int i = 0; i < _nss; ++i)
      _tau_b_tmp_old[i] = _tau_b_old[_qp][i];
  
  }
  else
  {
    if (_first_step_iter)
    {
      _gss_tmp = _gss_tmp_old = _gss_old[_qp];
      _accslip_tmp_old = _acc_slip_old[_qp];
	  
	  for (unsigned int i = 0; i < _nss; ++i)
	    _rho_for_tmp_old[i] = _rho_for_old[_qp][i];
  
      _rho_sub_tmp_old = _rho_sub_old[_qp];
	  
	  for (unsigned int i = 0; i < _nss; ++i)
        _tau_b_tmp_old[i] = _tau_b_old[_qp][i];
	
    }
    else
      _gss_tmp = _gss_tmp_old;
  }
}

// Assign updated temporary variables
// to the state variables
void
FiniteStrainCrystalPlasticityBackstress::postSolveStatevar()
{
  if (_max_substep_iter == 1) // No substepping
  {
    _gss[_qp] = _gss_tmp;
    _acc_slip[_qp] = _accslip_tmp;
	
	for (unsigned int i = 0; i < _nss; ++i)
	  _rho_for[_qp][i] = _rho_for_tmp[i];
  
    _rho_sub[_qp] = _rho_sub_tmp;
	
	for (unsigned int i = 0; i < _nss; ++i)
	  _tau_b[_qp][i] = _tau_b_tmp[i];
	
  }
  else
  {
    if (_last_step_iter)
    {
      _gss[_qp] = _gss_tmp;
      _acc_slip[_qp] = _accslip_tmp;
	  
	  for (unsigned int i = 0; i < _nss; ++i)
	    _rho_for[_qp][i] = _rho_for_tmp[i];
  
      _rho_sub[_qp] = _rho_sub_tmp;	

	  for (unsigned int i = 0; i < _nss; ++i)
	    _tau_b[_qp][i] = _tau_b_tmp[i];
  
    }
    else
    {
      _gss_tmp_old = _gss_tmp;
      _accslip_tmp_old = _accslip_tmp;

	  for (unsigned int i = 0; i < _nss; ++i)
	    _rho_for_tmp_old[i] = _rho_for_tmp[i];
  
      _rho_sub_tmp_old = _rho_sub_tmp;	 

      for (unsigned int i = 0; i < _nss; ++i)
	    _tau_b_tmp_old[i] = _tau_b_tmp[i];
	
    }
  }
}

void
FiniteStrainCrystalPlasticityBackstress::calcResidual( RankTwoTensor &resid )
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

  getSlipIncrements(); // Calculate dslip,dslipdtau, includes backstress

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
                      * ((1.0/2.0) * dCTE_dT * (temp-reference_temperature)*(temp-reference_temperature)
                      + thermal_expansion * (temp - reference_temperature))) - 1.0) * iden;
  pk2_new = _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);
  
  resid = _pk2_tmp - pk2_new;
  
  _lattice_strain[_qp] = _crysrot[_qp].transpose() * ee * _crysrot[_qp];
}

// Calculate slip increment,dslipdtau
// Critical resolved shear stress decreases exponentially with temperature
// A + B exp(- C * (T - 303.0))
// see equation (5) in:
// Nicolò Grilli , Alan C.F. Cocks , Edmund Tarleton
// Crystal plasticity finite element modelling of coarse-grained alpha-uranium
// Computational Materials Science 171 (2020) 109276
// Plus Armstrong-Frederick backstress term:
// Armstrong, P.J., Frederick, C.O., 1966. 
// A Mathematical Representation of the Multiaxial Bauschinger Effect, 
// G.E.G.B. Report RD/B/N. Central Electricity Generating Board.
void
FiniteStrainCrystalPlasticityBackstress::getSlipIncrements()
{
  Real temp = _temp[_qp];
  Real tau_b; // Backstress: temporary variable
  
  // Critical resolved shear stress in the input file
  // refers always to room temperature
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _gssT[i] = ( _dCRSS_dT_A + _dCRSS_dT_B * std::exp(- _dCRSS_dT_C * (temp - _reference_temperature))) * 
	           _gss_tmp[i];
  }
  
  for (unsigned int i = 0; i < _nss; ++i)
  {
	tau_b = _tau_b[_qp][i];
	
    _slip_incr(i) = _a0(i) * 
	                std::pow(std::abs((_tau(i) - tau_b) / _gssT[i]), 1.0 / _xm(i)) *
                    std::copysign(1.0, (_tau(i) - tau_b)) * _dt;
					
    if (std::abs(_slip_incr(i)) > _slip_incr_tol)
    {
      //_err_tol = true;
      //mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
	  
	  _slip_incr(i) = _slip_incr_tol * std::copysign(1.0, _tau(i));

    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
  {
	tau_b = _tau_b[_qp][i];
	  
    _dslipdtau(i) = _a0(i) / _xm(i) *
                    std::pow(std::abs((_tau(i) - tau_b) / _gssT[i]), 1.0 / _xm(i) - 1.0) / _gssT[i] *
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

// Calculate slip system resistance (CRSS)
// based on Taylor hardening model
// see equation (5) in:
// Nicolo Grilli , Alan C.F. Cocks , Edmund Tarleton
// Crystal plasticity finite element modelling of coarse-grained alpha-uranium
// Computational Materials Science 171 (2020) 109276
void
FiniteStrainCrystalPlasticityBackstress::updateGss()
{
  Real qab; // temporary variable to calculate Taylor hardening
  
  std::vector<Real> rho_for(_nss); // forest dislocation density
  Real rho_sub; // substructure dislocation density

  for (unsigned int i = 0; i < _nss; ++i)
    rho_for[i] = _rho_for[_qp][i]; // assign forest dislocation density
  
  rho_sub = _rho_sub[_qp]; // assign substructure dislocation density

  // store accumulated slip
  // not strictly necessary for the constitutive model
  _accslip_tmp = _accslip_tmp_old;
  for (unsigned int i = 0; i < _nss; ++i)
    _accslip_tmp += std::abs(_slip_incr(i));

  updateDisloDensity();
  
  updateBackstress();
    
  // qab must get the value in equation (5)
  // without the exponential of the temperature
  // because getSlipIncrements is taking care of that
  
  for (unsigned int i = 0; i < _nss; ++i)
  {
	qab = 0.0; // temporary variable
	qab += _tau0;
	qab += 0.9 * _burgers_vector * _projected_mu * std::sqrt(rho_for[i]);
	
	if(_rho_sub_flag) { // Model with rho_sub
	
	  qab -= 0.086 * _burgers_vector * _projected_mu * std::sqrt(rho_sub) 
	               * std::log(_burgers_vector * std::sqrt(rho_sub));
				   
	}
	
    _gss_tmp[i] = qab;
  }
}

// Calculate slip system resistance (CRSS)
// based on Taylor hardening model
// see equations (6) and (7) in:
// Nicolò Grilli , Alan C.F. Cocks , Edmund Tarleton
// Crystal plasticity finite element modelling of coarse-grained alpha-uranium
// Computational Materials Science 171 (2020) 109276
void
FiniteStrainCrystalPlasticityBackstress::updateDisloDensity()
{ 
  Real temp = _temp[_qp];
  Real da_tmp; // temporary variable to store the temperature dependent annihilation length	
  Real drho_for_tmp; // temporary variable to store the forest dislocation density increment
  Real drho_sub_tmp; // temporary variable to store the substructure dislocation density increment
  Real sum_rho_for; // sum of the forest dislocation density on all slip systems
  Real sum_slip_incr; // sum of the plastic slip increments on all slip systems
	
  for (unsigned int i = 0; i < _nss; ++i)	
    _rho_for_tmp[i] = _rho_for_tmp_old[i];

  _rho_sub_tmp = _rho_sub_tmp_old;
  
  sum_rho_for = 0.0;
  sum_slip_incr = 0.0;
  for (unsigned int i = 0; i < _nss; ++i) {
    sum_rho_for += _rho_for_tmp[i];
    sum_slip_incr += std::abs(_slip_incr(i)); // take abs because _slip_incr(i) is signed 
  }
  
  // equation (8) in the paper to define temperature dependent annihilation length 
  // note that in this implementation the annihilation length does not depend
  // on the specific slip system, it is realistic in FCC crystals  
  da_tmp = (_k_bol * temp) / (_drag_stress * std::pow(_burgers_vector,3.0));		 
  da_tmp *= _log_strain_rate_ratio; 
  da_tmp += 1.0;
  da_tmp *= _da0;
  
  // equation (6) in the paper
  // _dt is the substep of the crystal plasticity solver
  // rate equation must be multiplied by _dt to obtain time integration
  // otherwise use _slip_incr(i), which is the plastic strain rate
  // already multiplied by _dt  
  for (unsigned int i = 0; i < _nss; ++i) {
 
	drho_for_tmp = std::sqrt(_rho_for_tmp[i]) - da_tmp * _rho_for_tmp[i];
	drho_for_tmp *= _ka;
	
	// multiply by absolute value of _slip_incr(i)
	// because dislocation density growth is independent of the slip direction
	drho_for_tmp *= std::abs(_slip_incr(i));

	_rho_for_tmp[i] += drho_for_tmp;
  }    
  
  // equation (7) in the paper, but in this implementation not only the
  // first slip system contributes to the substructure dislocation density increase
  // but all the slip systems do, therefore multiply by sum_rho_for and
  if (_rho_sub_flag) { // Model with rho_sub
  
    drho_sub_tmp = 1800.0 * _ka * _burgers_vector * da_tmp * sum_rho_for * std::sqrt(_rho_sub_tmp);
    drho_sub_tmp *= sum_slip_incr;
	
  } else { // Model without rho_sub
	  
    drho_sub_tmp = 0.0;
	  
  }
  
  _rho_sub_tmp += drho_sub_tmp;
	
}

// Update of the Armstrong-Frederick backstress term:
// Armstrong, P.J., Frederick, C.O., 1966. 
// A Mathematical Representation of the Multiaxial Bauschinger Effect, 
// G.E.G.B. Report RD/B/N. Central Electricity Generating Board
// Backstress term can be changed by changing this function only
void
FiniteStrainCrystalPlasticityBackstress::updateBackstress()
{
  Real dtau_b_tmp; // increment of backstress, temporary variable
	
  for (unsigned int i = 0; i < _nss; ++i)
    _tau_b_tmp[i] = _tau_b_tmp_old[i];

  // _dt here is the substep
  for (unsigned int i = 0; i < _nss; ++i) {
	  
	dtau_b_tmp = 0.0;
	dtau_b_tmp = _c_backstress * _slip_incr(i);
	dtau_b_tmp -= _d_backstress * _tau_b_tmp[i] * std::abs(_slip_incr(i));
    
    _tau_b_tmp[i] += dtau_b_tmp;
  }
  
}

