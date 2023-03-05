// Nicolò Grilli
// University of Bristol
// 25 Gennaio 2023

#include "CrystalPlasticityDislocationTransport.h"
#include "libmesh/int_range.h"
#include <cmath>
#include "Function.h"

registerMooseObject("TensorMechanicsApp", CrystalPlasticityDislocationTransport);

InputParameters
CrystalPlasticityDislocationTransport::validParams()
{
  InputParameters params = CrystalPlasticityDislocationUpdateBase::validParams();
  params.addClassDescription("Dislocation based model for crystal plasticity with Orowan's law "
                             "using the stress update code. "
                             "Plastic slip is calculated using the CDD model variables, "
                             "therefore no update of material properties. ");
  params.addCoupledVar("temperature", 303.0,"Temperature, initialize at room temperature");
  params.addCoupledVar("rho_t_vector","Total dislocation density vector: each component is the total dislocation density on each slip system");
  params.addCoupledVar("rho_edge_vector","Total edge GND density vector: each component is the edge GND density on each slip system");
  params.addCoupledVar("rho_screw_vector","Total screw GND density vector: each component is the screw GND density on each slip system");
  params.addCoupledVar("q_t_vector","Curvature density vector: each component is the curvature on each slip system");
  params.addParam<Real>("burgers_vector_mag",0.000256,"Magnitude of the Burgers vector");
  params.addParam<Real>("shear_modulus",86000.0,"Shear modulus in Taylor hardening law G");
  params.addParam<Real>("alpha_0",0.4,"Prefactor of Taylor hardening law, alpha");
  params.addParam<Real>("r", 1.4, "Latent hardening coefficient");
  params.addParam<Real>("tau_c_0", 0.0, "Peierls stress");
  params.addParam<Real>("dislo_mobility",0.0,"Dislocation mobility");
  params.addParam<Real>("reduced_mobility",0.0,"Ratio between mobility above vmax and mobility");
  params.addParam<Real>("dislo_max_velocity",1000.0,"Maximum dislocation velocity (phonon drag)");
  params.addParam<Real>("bowout_coef",0.0,"bow-out coefficient: alpha in 4.30 of Hull-Bacon book");
  params.addParam<Real>("bowout_rho_threshold",0.2,"dislo density threshold to apply bow-out");
  params.addParam<Real>("rho_v_thres",0.001,"Dislo density threshold below which velocity goes to zero");
  params.addParam<bool>("rho_v_thres_flag",false,"Flag to determine whether to apply the previous threshold");
  params.addParam<MaterialPropertyName>(
      "total_twin_volume_fraction",
      "Total twin volume fraction, if twinning is considered in the simulation");
  params.addParam<UserObjectName>("read_initial_gnd_density",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element value "
                                  "of the initial GND density");
  params.addParam<Real>("reference_temperature",303.0,"reference temperature for thermal expansion");
  params.addParam<Real>("dCRSS_dT_A",1.0,"A coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 303.0))");
  params.addParam<Real>("dCRSS_dT_B",0.0,"B coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 303.0))");
  params.addParam<Real>("dCRSS_dT_C",0.0,"C coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 303.0))");
  return params;
}

CrystalPlasticityDislocationTransport::CrystalPlasticityDislocationTransport(
    const InputParameters & parameters)
  : CrystalPlasticityDislocationUpdateBase(parameters),
  
	// Temperature dependent properties
	_temperature(coupledValue("temperature")),

    // Dislocation model variables, not updated here
    _rho_t_vector(isParamValid("rho_t_vector") ? coupledVectorValue("rho_t_vector") : _vector_zero),	
    _rho_edge_vector(isParamValid("rho_edge_vector") ? coupledVectorValue("rho_edge_vector") : _vector_zero),
    _rho_screw_vector(isParamValid("rho_screw_vector") ? coupledVectorValue("rho_screw_vector") : _vector_zero),	
    _q_t_vector(isParamValid("q_t_vector") ? coupledVectorValue("q_t_vector") : _vector_zero),
  
    // Constitutive model parameters
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")),
	_shear_modulus(getParam<Real>("shear_modulus")),
	_alpha_0(getParam<Real>("alpha_0")),
    _r(getParam<Real>("r")),
	_tau_c_0(getParam<Real>("tau_c_0")),
	
	// Dislocation mobility parameters
	_dislo_mobility(getParam<Real>("dislo_mobility")),
	_reduced_mobility(getParam<Real>("reduced_mobility")),
    _dislo_max_velocity(getParam<Real>("dislo_max_velocity")), // Maximum dislocation velocity (phonon drag)
	_bowout_coef(getParam<Real>("bowout_coef")),
	_bowout_rho_threshold(getParam<Real>("bowout_rho_threshold")),
    _rho_v_thres(getParam<Real>("rho_v_thres")), // Dislo density threshold below which velocity goes to zero
    _rho_v_thres_flag(getParam<bool>("rho_v_thres_flag")), // Flag to determine whether to apply the previous threshold

    // Twinning contributions, if used
    _include_twinning_in_Lp(parameters.isParamValid("total_twin_volume_fraction")),
    _twin_volume_fraction_total(_include_twinning_in_Lp
                                    ? &getMaterialPropertyOld<Real>("total_twin_volume_fraction")
                                    : nullptr),
	
    // Temperature dependence of CRSS parameters
    _reference_temperature(getParam<Real>("reference_temperature")),
    _dCRSS_dT_A(getParam<Real>("dCRSS_dT_A")),
	_dCRSS_dT_B(getParam<Real>("dCRSS_dT_B")),
	_dCRSS_dT_C(getParam<Real>("dCRSS_dT_C")),
	
	// Dislocation velocity and its derivatives
	_dislo_velocity(declareProperty<std::vector<Real>>("dislo_velocity")), 
    _ddislo_velocity_dtau(declareProperty<std::vector<Real>>("ddislo_velocity_dtau")),

    // store edge and screw slip directions to calculate directional derivatives
    // of the plastic slip rate	
    _edge_slip_direction(declareProperty<std::vector<Real>>("edge_slip_direction")),
	_screw_slip_direction(declareProperty<std::vector<Real>>("screw_slip_direction"))
{
}

void
CrystalPlasticityDislocationTransport::initQpStatefulProperties()
{
  // Slip resistance is resized here
  CrystalPlasticityDislocationUpdateBase::initQpStatefulProperties();
  
  // Total dislocation density
  Real TotalRho = 0.0;
  
  // Temperature dependence of the CRSS
  Real temperature_dependence;
 
  // Critical resolved shear stress decreases exponentially with temperature
  // A + B exp(- C * (T - T0)) 
  temperature_dependence = ( _dCRSS_dT_A + _dCRSS_dT_B 
                         * std::exp(- _dCRSS_dT_C * (_temperature[_qp] - _reference_temperature)));
  
  // Initialize value of the slip resistance
  // as a function of the dislocation density
  for (const auto i : make_range(_number_slip_systems))
  {
    TotalRho += _rho_t_vector[_qp](i);
    
    // Add Peierls stress to each slip system
    _slip_resistance[_qp][i] = _tau_c_0;
  }
  
  if (TotalRho >= _bowout_rho_threshold) { // avoid that bow-out term becomes too large
	  
    for (const auto i : make_range(_number_slip_systems))
    {
	  // Taylor hardening + bow-out term
	  // See Hull, Bacon, Dislocations book equation 4.30
      _slip_resistance[_qp][i] += _alpha_0 * _shear_modulus * _burgers_vector_mag 
	    * (std::sqrt(TotalRho) + _bowout_coef * (_q_t_vector[_qp](i) / TotalRho));
	}
  
  } else if (TotalRho >= 0.0) {
	  
    for (const auto i : make_range(_number_slip_systems))
    {	  
      _slip_resistance[_qp][i] += _alpha_0 * _shear_modulus * _burgers_vector_mag * std::sqrt(TotalRho);	  
    }
  
  } else {
	  
    for (const auto i : make_range(_number_slip_systems))
    {
	  _slip_resistance[_qp][i] += 0.0;
    }
    
  }
  
  // Temperature dependence prefactor
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance[_qp][i] *= temperature_dependence;
  }

  // initialize slip increment
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_increment[_qp][i] = 0.0;
  }
  
  // Initialize vectors size here because they are used by AuxKernels
  // that are called just after initialization  
  _edge_slip_direction[_qp].resize(LIBMESH_DIM * _number_slip_systems);
  _screw_slip_direction[_qp].resize(LIBMESH_DIM * _number_slip_systems);

}

// Calculate Schmid tensor and
// store edge and screw slip directions to calculate directional derivatives
// of the plastic slip rate
void
CrystalPlasticityDislocationTransport::calculateSchmidTensor(
    const unsigned int & number_slip_systems,
    const std::vector<RealVectorValue> & plane_normal_vector,
    const std::vector<RealVectorValue> & direction_vector,
    std::vector<RankTwoTensor> & schmid_tensor,
    const RankTwoTensor & crysrot)
{
  std::vector<RealVectorValue> local_direction_vector, local_plane_normal;
  local_direction_vector.resize(number_slip_systems);
  local_plane_normal.resize(number_slip_systems);
  
  // Temporary directions and normals to calculate
  // screw dislocation slip direction
  RealVectorValue temp_mo;
  RealVectorValue temp_no;
  RealVectorValue temp_screw_mo;

  // Update slip direction and normal with crystal orientation
  for (const auto i : make_range(_number_slip_systems))
  {
    local_direction_vector[i].zero();
    local_plane_normal[i].zero();

    for (const auto j : make_range(LIBMESH_DIM))
      for (const auto k : make_range(LIBMESH_DIM))
      {
        local_direction_vector[i](j) =
            local_direction_vector[i](j) + crysrot(j, k) * direction_vector[i](k);

        local_plane_normal[i](j) =
            local_plane_normal[i](j) + crysrot(j, k) * plane_normal_vector[i](k);
      }

    // Calculate Schmid tensor
    for (const auto j : make_range(LIBMESH_DIM))
      for (const auto k : make_range(LIBMESH_DIM))
      {
        schmid_tensor[i](j, k) = local_direction_vector[i](j) * local_plane_normal[i](k);
      }
  }
  
  // Calculate and store edge and screw slip directions are also assigned
  _edge_slip_direction[_qp].resize(LIBMESH_DIM * _number_slip_systems);
  _screw_slip_direction[_qp].resize(LIBMESH_DIM * _number_slip_systems);
  
  for (const auto i : make_range(_number_slip_systems)) {
	for (const auto j : make_range(LIBMESH_DIM)) {
	  _edge_slip_direction[_qp][i * LIBMESH_DIM + j] = local_direction_vector[i](j);
	} 
  }
  
  for (const auto i : make_range(_number_slip_systems)) {
    for (const auto j : make_range(LIBMESH_DIM)) {
	  // assign temporary slip direction and normal for this slip system
      temp_mo(j) = local_direction_vector[i](j);
	  temp_no(j) = local_plane_normal[i](j);
    }  
	
	// calculate screw slip direction for this slip system
	// and store it in the screw slip direction vector
	temp_screw_mo = temp_mo.cross(temp_no);
	
	for (const auto j : make_range(LIBMESH_DIM)) {
	  _screw_slip_direction[_qp][i * LIBMESH_DIM + j] = temp_screw_mo(j);
	}
  }
  
}

// No costitutive variables here because dislocation density
// is a FE problem variable
// This is called in ComputeCrystalPlasticityStressDamage,
// therefore it must be kept here even if dummy.
void
CrystalPlasticityDislocationTransport::setInitialConstitutiveVariableValues()
{
  // Initialize state variables with the value at the previous time step
}

// This is called in ComputeCrystalPlasticityStressDamage,
// therefore it must be kept here even if dummy.
void
CrystalPlasticityDislocationTransport::setSubstepConstitutiveVariableValues()
{
  // Inialize state variable of the next substep
  // with the value at the previous substep
}

// Slip resistance can be calculated from dislocation density here only
// because it is the first method in which it is used,
// while calculateConstitutiveSlipDerivative is called afterwards.
// No backstress is included.
bool
CrystalPlasticityDislocationTransport::calculateSlipRate()
{
  calculateSlipResistance();
  
  // calculate dislocation velocity
  // and store it for advection kernel
  // necessary to call it here because slip rate
  // depends on dislocation velocity
  getDisloVelocity();

  // _slip_increment is the strain rate
  for (const auto i : make_range(_number_slip_systems))
  {
    if (_rho_t_vector[_qp](i) > 0.0) {
		
      _slip_increment[_qp][i] = _rho_t_vector[_qp](i) *
        std::abs(_dislo_velocity[_qp][i]) * _burgers_vector_mag *
        std::copysign(1.0, _tau[_qp][i]);
		
	} else {
		
      _slip_increment[_qp][i] = 0.0;
		
	}

    if (std::abs(_slip_increment[_qp][i]) * _substep_dt > _slip_incr_tol)
    {
      if (_print_convergence_message)
        mooseWarning("Maximum allowable slip increment exceeded ",
                     std::abs(_slip_increment[_qp][i]) * _substep_dt);

      return false;
    } 
  }
  
  return true;
}

// Slip resistance based on Taylor hardening
// and dislocation bow-out stress
void
CrystalPlasticityDislocationTransport::calculateSlipResistance()
{	
  // Total dislocation density
  Real TotalRho = 0.0;
  
  // Temperature dependence of the CRSS
  Real temperature_dependence;
  
  // Critical resolved shear stress decreases exponentially with temperature
  // A + B exp(- C * (T - T0)) 
  temperature_dependence = ( _dCRSS_dT_A + _dCRSS_dT_B 
                         * std::exp(- _dCRSS_dT_C * (_temperature[_qp] - _reference_temperature)));

  // Calculate value of the slip resistance
  // as a function of the dislocation density
  for (const auto i : make_range(_number_slip_systems))
  {
    TotalRho += _rho_t_vector[_qp](i);
    
    // Add Peierls stress to each slip system
    _slip_resistance[_qp][i] = _tau_c_0;
  }
  
  if (TotalRho >= _bowout_rho_threshold) { // avoid that bow-out term becomes too large
	  
    for (const auto i : make_range(_number_slip_systems))
    {
	  // Taylor hardening + bow-out term
	  // See Hull, Bacon, Dislocations book equation 4.30
      _slip_resistance[_qp][i] += _alpha_0 * _shear_modulus * _burgers_vector_mag 
	    * (std::sqrt(TotalRho) + _bowout_coef * (_q_t_vector[_qp](i) / TotalRho));
	}
  
  } else if (TotalRho >= 0.0) {
	  
    for (const auto i : make_range(_number_slip_systems))
    {	  
      _slip_resistance[_qp][i] += _alpha_0 * _shear_modulus * _burgers_vector_mag * std::sqrt(TotalRho);	  
    }
  
  } else {
	  
    for (const auto i : make_range(_number_slip_systems))
    {
	  _slip_resistance[_qp][i] += 0.0;
    }
    
  }
  
  // Temperature dependence prefactor
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance[_qp][i] *= temperature_dependence;
  }	

}

// Calculate dislocation velocity (edge and screw) as a function
// of the resolved shear stress and its derivative
void
CrystalPlasticityDislocationTransport::getDisloVelocity()
{		
  // resolved shear stress at max velocity _dislo_max_velocity	
  Real tau0; 

  _dislo_velocity[_qp].resize(_number_slip_systems);
  _ddislo_velocity_dtau[_qp].resize(_number_slip_systems);
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i) {
    _dislo_velocity[_qp][i] = 0.0;
    _ddislo_velocity_dtau[_qp][i] = 0.0;
  }
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
	tau0 = 0.0;
		
	if (_dislo_mobility > 0.0) {
	  tau0 = _dislo_max_velocity / _dislo_mobility; // temporary variable for this slip system
	  tau0 += _slip_resistance[_qp][i];		
	}
	
	if (std::abs(_tau[_qp][i]) > tau0) { // Case above _dislo_max_velocity: use reduced mobility
		
	  _dislo_velocity[_qp][i] = (_dislo_max_velocity + _reduced_mobility * (std::abs(_tau[_qp][i]) - tau0))
	                          * std::copysign(1.0, _tau[_qp][i]);
	  
	  // Derivative is always positive
	  _ddislo_velocity_dtau[_qp][i] = _reduced_mobility;
	  
	  if (_rho_v_thres_flag) { // Case with density below threshold

        if (_rho_t_vector[_qp](i) < _rho_v_thres) { // rescale dislocation velocity and derivative by a factor

          _dislo_velocity[_qp][i] *= (_rho_t_vector[_qp](i) / _rho_v_thres);
		  _ddislo_velocity_dtau[_qp][i] *= (_rho_t_vector[_qp](i) / _rho_v_thres);
		  
        }
		  
      }
		 
	} else if (std::abs(_tau[_qp][i]) > _slip_resistance[_qp][i]) { // Case below _dislo_max_velocity
		
      _dislo_velocity[_qp][i] = _dislo_mobility * (std::abs(_tau[_qp][i]) - _slip_resistance[_qp][i])
	                          * std::copysign(1.0, _tau[_qp][i]);

	  // Derivative is always positive
	  _ddislo_velocity_dtau[_qp][i] = _dislo_mobility;
		
	  if (_rho_v_thres_flag) { // Case with density below threshold

        if (_rho_t_vector[_qp](i) < _rho_v_thres) { // rescale dislocation velocity and derivative by a factor

          _dislo_velocity[_qp][i] *= (_rho_t_vector[_qp](i) / _rho_v_thres);
		  _ddislo_velocity_dtau[_qp][i] *= (_rho_t_vector[_qp](i) / _rho_v_thres);
		  
        }
		  
      }		
	  
	} else { // Case below critical resolved shear stress
		
	  _dislo_velocity[_qp][i] = 0.0;
	  _ddislo_velocity_dtau[_qp][i] = 0.0;
	  
	}	

  } // end cycle over slip systems
}

void
CrystalPlasticityDislocationTransport::calculateEquivalentSlipIncrement(
    RankTwoTensor & equivalent_slip_increment)
{
  if (_include_twinning_in_Lp)
  {
    for (const auto i : make_range(_number_slip_systems))
      equivalent_slip_increment += (1.0 - (*_twin_volume_fraction_total)[_qp]) *
                                   _flow_direction[_qp][i] * _slip_increment[_qp][i] * _substep_dt;
  }
  else // if no twinning volume fraction material property supplied, use base class
    CrystalPlasticityDislocationUpdateBase::calculateEquivalentSlipIncrement(equivalent_slip_increment);
}

// Note that this is always called after calculateSlipRate
// because calculateSlipRate is called in calculateResidual
// while this is called in calculateJacobian
// therefore it is ok to call calculateSlipResistance
// only inside calculateSlipRate
void
CrystalPlasticityDislocationTransport::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{	

  for (const auto i : make_range(_number_slip_systems))
  {
    dslip_dtau[i] = _rho_t_vector[_qp](i) *
      _ddislo_velocity_dtau[_qp][i] * _burgers_vector_mag;	  
  }	

}

// This is called in ComputeCrystalPlasticityStressDamage,
// therefore it must be kept here even if dummy.
bool
CrystalPlasticityDislocationTransport::areConstitutiveStateVariablesConverged()
{
  return true;						  
}

// This is called in ComputeCrystalPlasticityStressDamage,
// therefore it must be kept here even if dummy.
void
CrystalPlasticityDislocationTransport::updateSubstepConstitutiveVariableValues()
{
  // Update temporary variables at the end of the substep
}

// This is called in ComputeCrystalPlasticityStressDamage,
// therefore it must be kept here even if dummy.
void
CrystalPlasticityDislocationTransport::cacheStateVariablesBeforeUpdate()
{
  // Update the *_before_update variable
  // with internal state variables
}

// Dislocation density is a FE problem variable
// and not an internal variable,
// therefore, there is nothing to evolve here.
// This is called in ComputeCrystalPlasticityStressDamage,
// therefore it must be kept here even if dummy.
void
CrystalPlasticityDislocationTransport::calculateStateVariableEvolutionRateComponent()
{
}

// This is called in ComputeCrystalPlasticityStressDamage,
// therefore it must be kept here even if dummy.
bool
CrystalPlasticityDislocationTransport::updateStateVariables()
{
  // update internal state variables using _previous_substep_* variables
  // and perform the check to see if the state variable should be updated  
  return true;
}
