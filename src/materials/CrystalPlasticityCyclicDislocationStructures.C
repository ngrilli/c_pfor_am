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
  params.addParam<Real>("ao", 0.001, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");
  params.addParam<bool>("cap_slip_increment", false, "Cap the absolute value of the slip increment "
                                                     "in one time step to _slip_incr_tol. ");
  params.addParam<Real>("burgers_vector_mag",0.000256,"Magnitude of the Burgers vector");
  params.addParam<Real>("shear_modulus",86000.0,"Shear modulus in Taylor hardening law G");
  params.addParam<Real>("alpha_0",0.3,"Prefactor of Taylor hardening law, alpha");
  params.addParam<Real>("r", 1.4, "Latent hardening coefficient");
  params.addParam<Real>("tau_c_0", 0.112, "Peierls stress");
  params.addParam<Real>("k_0",100.0,"Coefficient K in SSD evolution, representing accumulation rate");
  params.addParam<Real>("y_c",0.0026,"Critical annihilation diameter");
  params.addParam<Real>("h",0.0,"Direct hardening coefficient for backstress");
  params.addParam<Real>("h_D",0.0,"Dynamic recovery coefficient for backstress");
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
	_alpha_0(getParam<Real>("alpha_0")),
    _r(getParam<Real>("r")),
	_tau_c_0(getParam<Real>("tau_c_0")),
	_k_0(getParam<Real>("k_0")),
	_y_c(getParam<Real>("y_c")),
	
	// Backstress parameters
	_h(getParam<Real>("h")),
	_h_D(getParam<Real>("h_D")),
	
	// Initial values of the state variables
    _init_rho_c(getParam<Real>("init_rho_c")),
	
	// State variables of the dislocation model
    _rho_c(declareProperty<std::vector<Real>>("rho_c")),
    _rho_c_old(getMaterialPropertyOld<std::vector<Real>>("rho_c")),
    
    // Walls fraction
    _f_w(declareProperty<Real>("f_w")),
    _f_w_old(declareProperty<Real>("f_w")),
    
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
  
  Real taylor_hardening;

  // Initialize the dislocation density size
  _rho_c[_qp].resize(_number_slip_systems);
  
  // Initialize the backstress size
  _backstress[_qp].resize(_number_slip_systems);
  
  // Initialize walls fraction
  _f_w[_qp] = 0.0;
  
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
    // Add Peierls stress
    _slip_resistance[_qp][i] = _tau_c_0;

    taylor_hardening = 0.0;
	  
    for (const auto j : make_range(_number_slip_systems))
    {
      // Determine slip planes
      unsigned int iplane, jplane;
      iplane = i / 3;
      jplane = j / 3;

      if (iplane == jplane) { // self vs. latent hardening
	  
	    // q_{ab} = 1.0 for self hardening
	    taylor_hardening += _rho_c[_qp][j]; 
		  
	  } else { // latent hardening
	  
	    taylor_hardening += (_r * _rho_c[_qp][j]);	  
		  
	  }
    }
	
	_slip_resistance[_qp][i] += (_alpha_0 * _shear_modulus * _burgers_vector_mag
	                          * std::sqrt(taylor_hardening));
	
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

  _f_w[_qp] = _f_w_old[_qp];

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
  Real taylor_hardening;

  for (const auto i : make_range(_number_slip_systems))
  {
    // Add Peierls stress
    _slip_resistance[_qp][i] = _tau_c_0;

    taylor_hardening = 0.0;
	  
    for (const auto j : make_range(_number_slip_systems))
    {
      // Determine slip planes
      unsigned int iplane, jplane;
      iplane = i / 3;
      jplane = j / 3;

      if (iplane == jplane) { // self vs. latent hardening
	  
	    // q_{ab} = 1.0 for self hardening
	    taylor_hardening += _rho_c[_qp][j]; 
		  
	  } else { // latent hardening
	  
	    taylor_hardening += (_r * _rho_c[_qp][j]);

	  }
    }
	
	_slip_resistance[_qp][i] += (_alpha_0 * _shear_modulus * _burgers_vector_mag
	                          * std::sqrt(taylor_hardening));
	
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
  Real rho_sum;

  // SSD dislocation density increment
  for (const auto i : make_range(_number_slip_systems))
  {
    
    rho_sum = _rho_c[_qp][i];

    // Multiplication and annihilation
	// note that _slip_increment here is the rate
	// and the rate equation gets multiplied by time step in updateStateVariables
    _rho_c_increment[i] = _k_0 * sqrt(rho_sum) - 2 * _y_c * _rho_c[_qp][i];
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
    _backstress_increment[i] = _h * _slip_increment[_qp][i];
    _backstress_increment[i] -= _h_D * _backstress[_qp][i] * std::abs(_slip_increment[_qp][i]);  
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
