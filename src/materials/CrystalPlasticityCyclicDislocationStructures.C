// Marco Pelegatti
// Università di Udine
// Nicolò Grilli
// Università di Bristol
// 23 Marzo 2024

#include "CrystalPlasticityCyclicDislocationStructures.h"
#include "libmesh/int_range.h"
#include <cmath>
#include "Function.h"

registerMooseObject("c_pfor_amApp", CrystalPlasticityCyclicDislocationStructures);

InputParameters
CrystalPlasticityCyclicDislocationStructures::validParams()
{
  InputParameters params = CrystalPlasticityDislocationUpdateBase::validParams();
  params.addClassDescription("Dislocation based model for crystal plasticity "
                             "using the stress update code. "
                             "Includes slip, creep and backstress. ");
  params.addParam<Real>("ao", 0.004, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");
  params.addParam<bool>("cap_slip_increment", false, "Cap the absolute value of the slip increment "
                                                     "in one time step to _slip_incr_tol");
  params.addParam<Real>("burgers_vector_mag", 0.000256, "Magnitude of the Burgers vector");
  params.addParam<Real>("shear_modulus", 76800.0, "Shear modulus in Taylor hardening law G");
  params.addParam<Real>("A_self", 0.122, "Self hardening coefficient");
  params.addParam<Real>("nu", 0.33, "Poisson's ratio for backstress calculation by Eshelby's inclusion");
  params.addParam<Real>("K_struct",3.0,"Constant of similitude for dislocation substructure");
  params.addParam<Real>("tau_0", 80.0, "Resolved shear stress at initial yield");
  params.addParam<Real>("k_0",100.0,"Coefficient K in SSD evolution, representing accumulation rate");
  params.addParam<Real>("y_c",0.013,"Critical annihilation diameter");
  params.addParam<Real>("f_0",0.5,"Initial dislocation walls volume fraction");
  params.addParam<Real>("f_inf",0.25,"Saturated dislocation walls volume fraction");
  params.addParam<Real>("k_f",0.5,"Rate constant for the evolution of the dislocation walls volume fraction");
  params.addParam<Real>("eta_0",50,"Initial Max/Min axis length ratio of the dislocation substructure");
  params.addParam<Real>("eta_inf",1,"Saturated Max/Min axis length ratio of the dislocation substructure");
  params.addParam<Real>("X_norm",400,"Normalization constant for the evolution of the dislocation substructure");
  params.addParam<Real>("init_d_struct",10.0,"Initial characteristic dislocation substructure length");
  params.addParam<Real>("init_rho_c",1.0,"Initial channel dislocation density");
  return params;
}

CrystalPlasticityCyclicDislocationStructures::CrystalPlasticityCyclicDislocationStructures(
    const InputParameters & parameters)
  : CrystalPlasticityDislocationUpdateBase(parameters),
  
    // Constitutive model parameters
    _ao(getParam<Real>("ao")),
    _xm(getParam<Real>("xm")),
    _cap_slip_increment(getParam<bool>("cap_slip_increment")),
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")),
	_shear_modulus(getParam<Real>("shear_modulus")),
	_A_self(getParam<Real>("A_self")),
	_nu(getParam<Real>("nu")),
	_K_struct(getParam<Real>("K_struct")),
	_tau_0(getParam<Real>("tau_0")),
	_k_0(getParam<Real>("k_0")),
	_y_c(getParam<Real>("y_c")),
	_f_0(getParam<Real>("f_0")),
	_f_inf(getParam<Real>("f_inf")),
	_k_f(getParam<Real>("k_f")),
	_eta_0(getParam<Real>("eta_0")),
	_eta_inf(getParam<Real>("eta_inf")),
	_X_norm(getParam<Real>("X_norm")),
	_init_d_struct(getParam<Real>("init_d_struct")),

	// Initial values of the state variables
    _init_rho_c(getParam<Real>("init_rho_c")),
	
	// State variables of the dislocation model
    _rho_c(declareProperty<std::vector<Real>>("rho_c")),
    _rho_c_old(getMaterialPropertyOld<std::vector<Real>>("rho_c")),
    
	// Cumulative effective plastic strain
	_epsilon_p_eff_cum(getMaterialProperty<Real>("epsilon_p_eff_cum")),    
    
	// Instantaneous plastic deformation tangent at the slip system level
	_dslip_dtau(declareProperty<std::vector<Real>>("dslip_dtau")),
	
    // Dislocation walls volume fraction
    _f_w(declareProperty<Real>("f_w")),

	// Max/min axis length ratio of the dislocation substructure
    _eta(declareProperty<Real>("eta")),
	
	// Characteristic dislocation substructure length
    _d_struct(declareProperty<Real>("d_struct")),
	
	// Mean glide distance for dislocations in the channel phase
    _l_c(declareProperty<Real>("l_c")),    
    
    // Backstress variable
    _backstress(declareProperty<std::vector<Real>>("backstress")),
    _backstress_old(getMaterialPropertyOld<std::vector<Real>>("backstress")),

    // increments of state variables
    _rho_c_increment(_number_slip_systems, 0.0),
    _backstress_increment(_number_slip_systems, 0.0),
	
	// resize local caching vectors used for substepping
	_previous_substep_rho_c(_number_slip_systems, 0.0),
	_previous_substep_backstress(_number_slip_systems, 0.0),
    _rho_c_before_update(_number_slip_systems, 0.0),
    _backstress_before_update(_number_slip_systems, 0.0)
{
}

void
CrystalPlasticityCyclicDislocationStructures::initQpStatefulProperties()
{
  // Slip resistance is resized here
  CrystalPlasticityDislocationUpdateBase::initQpStatefulProperties();

  // Initialize the dislocation density size
  _rho_c[_qp].resize(_number_slip_systems);
  
  // Initialize the backstress size
  _backstress[_qp].resize(_number_slip_systems);
  
  // Initialize instantaneous plastic deformation tangent at the slip system level
  _dslip_dtau[_qp].resize(_number_slip_systems);
  
  // Initialize walls fraction
  _f_w[_qp] = _f_0;
  
  // Initialize max/min axis length ratio of the dislocation substructure
  _eta[_qp] = _eta_0;
  
  // Initialize characteristic dislocation substructure length
  _d_struct[_qp] = _init_d_struct;
  
  // Initialize mean glide distance for dislocations in the channel phase
  _l_c[_qp] = _eta_0 * _init_d_struct;
  
  // Initialize dislocation densities and backstress
  for (const auto i : make_range(_number_slip_systems))
  {
    _rho_c[_qp][i] = _init_rho_c;
	
	_backstress[_qp][i] = 0.0;	
  }

  // Initialize value of the slip resistance
  // as a function of the dislocation density
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance[_qp][i] = _shear_modulus * _burgers_vector_mag
	                           * std::sqrt(_A_self * _rho_c[_qp][i]);
  }

  // initialize slip increment
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_increment[_qp][i] = 0.0;
  }
}

void
CrystalPlasticityCyclicDislocationStructures::setInitialConstitutiveVariableValues()
{
  // Initialize state variables with the value at the previous time step
  _rho_c[_qp] = _rho_c_old[_qp];
  _previous_substep_rho_c = _rho_c_old[_qp];
  _backstress[_qp] = _backstress_old[_qp];
  _previous_substep_backstress = _backstress_old[_qp];
}

void
CrystalPlasticityCyclicDislocationStructures::setSubstepConstitutiveVariableValues()
{
  // Inialize state variable of the next substep
  // with the value at the previous substep
  _rho_c[_qp] = _previous_substep_rho_c;
  _backstress[_qp] = _previous_substep_backstress;
}

// Slip resistance can be calculated from dislocation density here only
// because it is the first method in which it is used,
// while calculateConstitutiveSlipDerivative is called afterwards
bool
CrystalPlasticityCyclicDislocationStructures::calculateSlipRate()
{
  calculateSlipResistance();
  
  // Ratio between effective stress and CRSS
  // temporary variable for each slip system
  Real stress_ratio;
  
  // Difference between RSS and backstress
  // temporary variable for each slip system
  Real effective_stress;
  
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress = _tau[_qp][i] - _backstress[_qp][i];
    
    stress_ratio = std::abs(effective_stress / _slip_resistance[_qp][i]);
    
    _slip_increment[_qp][i] =
        _ao * std::pow(stress_ratio, 1.0 / _xm);
      
    if (effective_stress < 0.0)
      _slip_increment[_qp][i] *= -1.0;

    if (std::abs(_slip_increment[_qp][i]) * _substep_dt > _slip_incr_tol)
    {
      if (_cap_slip_increment) {
		  
        _slip_increment[_qp][i] = _slip_incr_tol * std::copysign(1.0, _slip_increment[_qp][i])
                                / _substep_dt;
                                
	  } else if (_print_convergence_message) {
		  
        mooseWarning("Maximum allowable slip increment exceeded ",
                     std::abs(_slip_increment[_qp][i]) * _substep_dt);

        return false;
	  }
    }
  }
  return true;
}

// Slip resistance based on Taylor hardening
void
CrystalPlasticityCyclicDislocationStructures::calculateSlipResistance()
{
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance[_qp][i] = _shear_modulus * _burgers_vector_mag
	                           * std::sqrt(_A_self * _rho_c[_qp][i]);
  }
}

// Note that this is always called after calculateSlipRate
// because calculateSlipRate is called in calculateResidual
// while this is called in calculateJacobian
// therefore it is ok to calculate calculateSlipRate
// only inside calculateSlipRate
void
CrystalPlasticityCyclicDislocationStructures::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  // Ratio between effective stress and CRSS
  // temporary variable for each slip system
  Real stress_ratio;
  
  // Difference between RSS and backstress
  // temporary variable for each slip system
  Real effective_stress;
	
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress = _tau[_qp][i] - _backstress[_qp][i];	  
	  
    if (MooseUtils::absoluteFuzzyEqual(effective_stress, 0.0)) {
		
      dslip_dtau[i] = 0.0;
      		
	} else {
		
	  stress_ratio = std::abs(effective_stress / _slip_resistance[_qp][i]);

      dslip_dtau[i] = _ao / _xm *
                      std::pow(stress_ratio, 1.0 / _xm - 1.0) /
                      _slip_resistance[_qp][i];		
	}
  }
}

void
CrystalPlasticityCyclicDislocationStructures::updateSubstepConstitutiveVariableValues()
{
  // Update temporary variables at the end of the substep
  _previous_substep_rho_c = _rho_c[_qp];
  _previous_substep_backstress = _backstress[_qp];
}

void
CrystalPlasticityCyclicDislocationStructures::cacheStateVariablesBeforeUpdate()
{
  _rho_c_before_update = _rho_c[_qp];
  _backstress_before_update = _backstress[_qp];
}

void
CrystalPlasticityCyclicDislocationStructures::calculateStateVariableEvolutionRateComponent()
{
  // channel dislocation density increment
  for (const auto i : make_range(_number_slip_systems))
  {
    // Multiplication and annihilation
	// note that _slip_increment here is the rate
	// and the rate equation gets multiplied by time step in updateStateVariables
    _rho_c_increment[i] = _k_0 * sqrt(_rho_c[_qp][i]) - 2 * _y_c * _rho_c[_qp][i];
    _rho_c_increment[i] *= std::abs(_slip_increment[_qp][i]) / _burgers_vector_mag;
  }
  
  // backstress increment
  ArmstrongFrederickBackstressUpdate();
}

// Armstrong-Frederick update of the backstress
void
CrystalPlasticityCyclicDislocationStructures::ArmstrongFrederickBackstressUpdate()
{
  for (const auto i : make_range(_number_slip_systems)) 
  {
    _backstress_increment[i] = _slip_increment[_qp][i];
    _backstress_increment[i] -= _backstress[_qp][i] * std::abs(_slip_increment[_qp][i]);  
  }
}

bool
CrystalPlasticityCyclicDislocationStructures::updateStateVariables()
{
  // Now perform the check to see if the slip system should be updated
  // SSD
  for (const auto i : make_range(_number_slip_systems))
  { 
    _rho_c_increment[i] *= _substep_dt;

    // force positive SSD density
    if (_previous_substep_rho_c[i] < _zero_tol && _rho_c_increment[i] < 0.0)
      _rho_c[_qp][i] = _previous_substep_rho_c[i];
    else
      _rho_c[_qp][i] = _previous_substep_rho_c[i] + _rho_c_increment[i];

    if (_rho_c[_qp][i] < 0.0)
      return false;
  }
  
  // Backstress: can be both positive or negative
  for (const auto i : make_range(_number_slip_systems))
  { 
    _backstress_increment[i] *= _substep_dt;
    _backstress[_qp][i] = _previous_substep_backstress[i] + _backstress_increment[i];
  }
  
  return true;
}
