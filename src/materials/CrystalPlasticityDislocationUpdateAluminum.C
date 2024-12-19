// Daijun Hu
// National University of Singapore
// 22 May 2024

#include "CrystalPlasticityDislocationUpdateAluminum.h"
#include "libmesh/int_range.h"
#include <cmath>
#include "Function.h"
#include "ComputeElasticityTensorCPGrain.h"

registerMooseObject("c_pfor_amApp", CrystalPlasticityDislocationUpdateAluminum);

InputParameters
CrystalPlasticityDislocationUpdateAluminum::validParams()
{
  InputParameters params = CrystalPlasticityDislocationUpdateBase::validParams();
  params.addClassDescription("Dislocation based model for crystal plasticity "
                             "using the stress update code. "
                             "Includes slip, creep and backstress. Modified for AM aluminum alloys "
                             "using a specific temperature dependence of the CRSS and strain rate sensitivity. ");
  params.addParam<Real>("ao", 0.001, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");  
  params.addParam<MaterialPropertyName>("xm_matprop",
    "Optional xm material property for exponent for slip rate. ");
  params.addParam<Real>("creep_ao", 0.0, "creep rate coefficient");
  params.addParam<Real>("creep_xm", 0.1, "exponent for creep rate");
  params.addParam<Real>("xm_max", 0.2, "upper limit of xm at semi-solid state");  
  params.addParam<Real>("xm_cali", 1.5, "coefficient to calibrate xm at initial state");  
  params.addParam<FunctionName>("creep_ao_function",
    "Optional function for creep prefactor. If provided, the creep prefactor can be set as a function of time. "
    "This is useful for an initial plastic deformation followed by creep load. ");
  params.addParam<Real>("burgers_vector_mag",0.000256,"Magnitude of the Burgers vector");
  params.addParam<Real>("shear_modulus",86000.0,"Shear modulus in Taylor hardening law G");
  params.addParam<Real>("alpha_0",0.3,"Prefactor of Taylor hardening law, alpha");
  params.addParam<Real>("r", 1.4, "Latent hardening coefficient");
  params.addParam<Real>("tau_c_0", 0.112, "Peierls stress");
  params.addParam<Real>("k_0",100.0,"Coefficient K in SSD evolution, representing accumulation rate");
  params.addParam<Real>("y_c",0.0026,"Critical annihilation diameter");
  params.addParam<Real>("h",0.0,"Direct hardening coefficient for backstress");
  params.addParam<Real>("h_D",0.0,"Dynamic recovery coefficient for backstress");
  params.addParam<Real>("rho_tol",1.0,"Tolerance on dislocation density update");
  params.addParam<Real>("init_rho_ssd",1.0,"Initial dislocation density");
  params.addParam<Real>("init_rho_gnd_edge",0.0,"Initial dislocation density");
  params.addParam<Real>("init_rho_gnd_screw",0.0,"Initial dislocation density");
  params.addParam<MaterialPropertyName>(
      "total_twin_volume_fraction",
      "Total twin volume fraction, if twinning is considered in the simulation");
  params.addParam<UserObjectName>("read_initial_gnd_density",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element value "
                                  "of the initial GND density");
  params.addCoupledVar("dslip_increment_dedge",0.0,"Directional derivative of the slip rate along the edge motion direction.");
  params.addCoupledVar("dslip_increment_dscrew",0.0,"Directional derivative of the slip rate along the screw motion direction.");
  params.addCoupledVar("temperature", 303.0,"Temperature, initialize at room temperature");
  params.addParam<Real>("reference_temperature",303.0,"reference temperature for thermal expansion");
  params.addParam<Real>("k_A",1.0,"k_A coefficient for the Boltzmann dependence of the critical "
                        "resolved shear stress with temperature. ");
  params.addParam<Real>("k_B",0.0,"k_B coefficient for the Boltzmann dependence of the critical "
                        "resolved shear stress with temperature. ");
  params.addParam<Real>("k_C",0.0,"k_C coefficient for the Boltzmann dependence of the critical "
                        "resolved shear stress with temperature. ");
  params.addParam<Real>("k_D",0.0,"k_D coefficient for the Boltzmann dependence of the critical "
                        "resolved shear stress with temperature. ");     
  params.addParam<Real>("melting_temperature_high", 925, "Melting temperature (liquidus) = zero stiffness.");  
  params.addParam<Real>("melting_temperature_low", 806, "Solidus = full stiffness.");                                          
  return params;
}

CrystalPlasticityDislocationUpdateAluminum::CrystalPlasticityDislocationUpdateAluminum(
    const InputParameters & parameters)
  : CrystalPlasticityDislocationUpdateBase(parameters),

    // Constitutive model parameters
    _ao(getParam<Real>("ao")),
    _xm(getParam<Real>("xm")),
    _include_xm_matprop(parameters.isParamValid("xm_matprop")),
    _xm_matprop(_include_xm_matprop
                ? &getMaterialProperty<Real>("xm_matprop")
                : nullptr),
    _creep_ao(getParam<Real>("creep_ao")),
    _creep_xm(getParam<Real>("creep_xm")),
    _creep_ao_function(this->isParamValid("creep_ao_function")
                       ? &this->getFunction("creep_ao_function")
                       : NULL),
    _xm_max(getParam<Real>("xm_max")),  
    _xm_cali(getParam<Real>("xm_cali")),      
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
	
	// Tolerance on dislocation density update
	_rho_tol(getParam<Real>("rho_tol")),	
	
	// Initial values of the state variables
    _init_rho_ssd(getParam<Real>("init_rho_ssd")),
    _init_rho_gnd_edge(getParam<Real>("init_rho_gnd_edge")),
    _init_rho_gnd_screw(getParam<Real>("init_rho_gnd_screw")),

    // Elasticity tensor used for strain rate sensitivity formula
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>("elasticity_tensor")),
	
	// State variables of the dislocation model
    _rho_ssd(declareProperty<std::vector<Real>>("rho_ssd")),
    _rho_ssd_old(getMaterialPropertyOld<std::vector<Real>>("rho_ssd")),
    _rho_gnd_edge(declareProperty<std::vector<Real>>("rho_gnd_edge")),
   	_rho_gnd_edge_old(getMaterialPropertyOld<std::vector<Real>>("rho_gnd_edge")),
  	_rho_gnd_screw(declareProperty<std::vector<Real>>("rho_gnd_screw")),
    _rho_gnd_screw_old(getMaterialPropertyOld<std::vector<Real>>("rho_gnd_screw")),
    
    // Backstress variable
    _backstress(declareProperty<std::vector<Real>>("backstress")),
    _backstress_old(getMaterialPropertyOld<std::vector<Real>>("backstress")),

    // increments of state variables
    _rho_ssd_increment(_number_slip_systems, 0.0),
    _rho_gnd_edge_increment(_number_slip_systems, 0.0),
    _rho_gnd_screw_increment(_number_slip_systems, 0.0),
    _backstress_increment(_number_slip_systems, 0.0),
	
	// resize local caching vectors used for substepping
    _previous_substep_rho_ssd(_number_slip_systems, 0.0),
	_previous_substep_rho_gnd_edge(_number_slip_systems, 0.0),
	_previous_substep_rho_gnd_screw(_number_slip_systems, 0.0),
	_previous_substep_backstress(_number_slip_systems, 0.0),
    _rho_ssd_before_update(_number_slip_systems, 0.0),
    _rho_gnd_edge_before_update(_number_slip_systems, 0.0),
    _rho_gnd_screw_before_update(_number_slip_systems, 0.0),  	
    _backstress_before_update(_number_slip_systems, 0.0),

    // Twinning contributions, if used
    _include_twinning_in_Lp(parameters.isParamValid("total_twin_volume_fraction")),
    _twin_volume_fraction_total(_include_twinning_in_Lp
                                    ? &getMaterialPropertyOld<Real>("total_twin_volume_fraction")
                                    : nullptr),

    // UserObject to read the initial GND density from file						
    _read_initial_gnd_density(isParamValid("read_initial_gnd_density")
                               ? &getUserObject<ElementPropertyReadFile>("read_initial_gnd_density")
                               : nullptr),
									
    // Directional derivatives of the slip rate
    _dslip_increment_dedge(coupledArrayValue("dslip_increment_dedge")), 
    _dslip_increment_dscrew(coupledArrayValue("dslip_increment_dscrew")),
	
	// Temperature dependent properties
	_temperature(coupledValue("temperature")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    // parameters for the Boltzmann function for the temperature dependence of the CRSS
    _k_A(getParam<Real>("k_A")),
	_k_B(getParam<Real>("k_B")),
	_k_C(getParam<Real>("k_C")),
    _k_D(getParam<Real>("k_D")),
    
    // store edge and screw slip directions to calculate directional derivatives
    // of the plastic slip rate	
    _edge_slip_direction(declareProperty<std::vector<Real>>("edge_slip_direction")),
	_screw_slip_direction(declareProperty<std::vector<Real>>("screw_slip_direction")),
  	_melting_temperature_high(getParam<Real>("melting_temperature_high")),
	_melting_temperature_low(getParam<Real>("melting_temperature_low"))
{
}

void
CrystalPlasticityDislocationUpdateAluminum::initQpStatefulProperties()
{
  // Slip resistance is resized here
  CrystalPlasticityDislocationUpdateBase::initQpStatefulProperties();
  
  Real taylor_hardening;
  
  // Temperature dependence of the CRSS
  Real temperature_dependence;
  Real temperature_Boltzmann; 

  // Initialize the dislocation density size
  _rho_ssd[_qp].resize(_number_slip_systems);
  _rho_gnd_edge[_qp].resize(_number_slip_systems);
  _rho_gnd_screw[_qp].resize(_number_slip_systems);
  
  // Initialize the backstress size
  _backstress[_qp].resize(_number_slip_systems);
  
  // Initialize dislocation densities and backstress
  for (const auto i : make_range(_number_slip_systems))
  {
    _rho_ssd[_qp][i] = _init_rho_ssd;
	
	if (_read_initial_gnd_density) { // Read initial GND density from file
	
    _rho_gnd_edge[_qp][i] = _read_initial_gnd_density->getData(_current_elem, i);
	_rho_gnd_screw[_qp][i] = _read_initial_gnd_density->getData(_current_elem, _number_slip_systems+i);
	
	} else { // Initialize uniform GND density
		
    _rho_gnd_edge[_qp][i] = _init_rho_gnd_edge;
    _rho_gnd_screw[_qp][i] = _init_rho_gnd_screw;
	
	}
	
	_backstress[_qp][i] = 0.0;
  }
 
  // Critical resolved shear stress decreases based on Boltzmann dependence with temperature
  temperature_Boltzmann = 1.0 + std::exp((_temperature[_qp] - _reference_temperature - _k_C) / _k_D);
  temperature_dependence = ((_k_A - _k_B) / temperature_Boltzmann) + _k_B;
  
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
	    taylor_hardening += (_rho_ssd[_qp][j] 
		          + std::abs(_rho_gnd_edge[_qp][j])
				  + std::abs(_rho_gnd_screw[_qp][j])); 
		  
	  } else { // latent hardening
	  
	    taylor_hardening += (_r * (_rho_ssd[_qp][j] 
		          + std::abs(_rho_gnd_edge[_qp][j])
				  + std::abs(_rho_gnd_screw[_qp][j])));	  
		  
	  }
    }
	
	_slip_resistance[_qp][i] += (_alpha_0 * _shear_modulus * _burgers_vector_mag
	                          * std::sqrt(taylor_hardening) * temperature_dependence);
	
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
CrystalPlasticityDislocationUpdateAluminum::calculateSchmidTensor(
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

void
CrystalPlasticityDislocationUpdateAluminum::setInitialConstitutiveVariableValues()
{
  // Initialize state variables with the value at the previous time step
  _rho_ssd[_qp] = _rho_ssd_old[_qp];
  _previous_substep_rho_ssd = _rho_ssd_old[_qp];
  _rho_gnd_edge[_qp] = _rho_gnd_edge_old[_qp];
  _previous_substep_rho_gnd_edge = _rho_gnd_edge_old[_qp];
  _rho_gnd_screw[_qp] = _rho_gnd_screw_old[_qp];
  _previous_substep_rho_gnd_screw = _rho_gnd_screw_old[_qp];
  _backstress[_qp] = _backstress_old[_qp];
  _previous_substep_backstress = _backstress_old[_qp];
}

void
CrystalPlasticityDislocationUpdateAluminum::setSubstepConstitutiveVariableValues()
{
  // Inialize state variable of the next substep
  // with the value at the previous substep
  _rho_ssd[_qp] = _previous_substep_rho_ssd;
  _rho_gnd_edge[_qp] = _previous_substep_rho_gnd_edge;
  _rho_gnd_screw[_qp] = _previous_substep_rho_gnd_screw;
  _backstress[_qp] = _previous_substep_backstress;
}

// Slip resistance can be calculated from dislocation density here only
// because it is the first method in which it is used,
// while calculateConstitutiveSlipDerivative is called afterwards
bool
CrystalPlasticityDislocationUpdateAluminum::calculateSlipRate()
{
  calculateSlipResistance();
  
  // Ratio between effective stress and CRSS
  // temporary variable for each slip system
  Real stress_ratio;
  
  // Difference between RSS and backstress
  // temporary variable for each slip system
  Real effective_stress;
  
  // Creep prefactor: if function is not given
  // the constant value is used
  Real creep_ao;
  
  // Strain rate sensitivity: if material property is not given
  // the temperature dependent value is used
  Real xm;
  Real mu_f;
  Real eta_f;
  Real k_b = 1.38e-11;
  Real xm_temp;
  
  if (_creep_ao_function) {
	  
    creep_ao = _creep_ao_function->value(_t, _q_point[_qp]);
	  
  } else {
	  
    creep_ao = _creep_ao;
	  
  }
   
  if (_include_xm_matprop) {
	  
    xm = (*_xm_matprop)[_qp];
	  
  } else {
	  
    mu_f = _elasticity_tensor[_qp](1, 2, 1, 2);
    eta_f = 0.5 * (_elasticity_tensor[_qp](0, 0, 0, 0) -  _elasticity_tensor[_qp](0, 0, 1, 1));
    xm_temp = 9 *_xm_cali * k_b *_temperature[_qp] / (_burgers_vector_mag * _burgers_vector_mag * _burgers_vector_mag * sqrt(mu_f * eta_f));
  
    if (_temperature[_qp] < _melting_temperature_high) {
    
      if (xm_temp <= _xm_max) {

        if (xm_temp >= _xm) {
          xm = xm_temp; 
        } else {
          xm = _xm;
        }
      } else {
        xm = _xm_max;
      }
    } else { // temperature above melting point
      xm = _xm;
    }
  }
  
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress = _tau[_qp][i] - _backstress[_qp][i];
    
    stress_ratio = std::abs(effective_stress / _slip_resistance[_qp][i]);
    
    _slip_increment[_qp][i] =
        _ao * std::pow(stress_ratio, 1.0 / xm)
      + creep_ao * std::pow(stress_ratio, 1.0 / _creep_xm);
      
    if (effective_stress < 0.0)
      _slip_increment[_qp][i] *= -1.0;

    if (std::abs(_slip_increment[_qp][i]) * _substep_dt > _slip_incr_tol)
    {

      _slip_increment[_qp][i] = _slip_incr_tol * std::copysign(1.0, _slip_increment[_qp][i]) / _substep_dt;

    }
  }
  
  return true;
}

// Slip resistance based on Taylor hardening
void
CrystalPlasticityDislocationUpdateAluminum::calculateSlipResistance()
{
  Real taylor_hardening;
  
  // Temperature dependence of the CRSS
  Real temperature_dependence;
  Real temperature_Boltzmann;  
  
  // Critical resolved shear stress decreases based on Boltzmann dependence with temperature
  temperature_Boltzmann = 1.0 + std::exp((_temperature[_qp] - _reference_temperature - _k_C) / _k_D);
  temperature_dependence = ((_k_A - _k_B) / temperature_Boltzmann) + _k_B;
	
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
	    taylor_hardening += (_rho_ssd[_qp][j] 
		          + std::abs(_rho_gnd_edge[_qp][j])
				  + std::abs(_rho_gnd_screw[_qp][j])); 
		  
	  } else { // latent hardening
	  
	    taylor_hardening += (_r * (_rho_ssd[_qp][j] 
		          + std::abs(_rho_gnd_edge[_qp][j])
				  + std::abs(_rho_gnd_screw[_qp][j])));
	  }
    }
	
	_slip_resistance[_qp][i] += (_alpha_0 * _shear_modulus * _burgers_vector_mag
	                          * std::sqrt(taylor_hardening) * temperature_dependence);
  }
}

void
CrystalPlasticityDislocationUpdateAluminum::calculateEquivalentSlipIncrement(
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
// therefore it is ok to calculate calculateSlipRate
// only inside calculateSlipRate
void
CrystalPlasticityDislocationUpdateAluminum::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  // Ratio between effective stress and CRSS
  // temporary variable for each slip system
  Real stress_ratio;
  
  // Difference between RSS and backstress
  // temporary variable for each slip system
  Real effective_stress;
  
  // Creep prefactor: if function is not given
  // the constant value is used
  Real creep_ao;
  
  // Strain rate sensitivity: if material property is not given
  // the temperature dependent value is used
  Real xm;
  Real mu_f;
  Real eta_f;
  Real k_b = 1.38e-11;
  Real xm_temp;
  
  if (_creep_ao_function) {
	  
    creep_ao = _creep_ao_function->value(_t, _q_point[_qp]);
	  
  } else {
	  
    creep_ao = _creep_ao;
	  
  }	

  if (_include_xm_matprop) {
	  
    xm = (*_xm_matprop)[_qp];	  
	  
  } else {
	  
    mu_f = _elasticity_tensor[_qp](1, 2, 1, 2);
    eta_f = 0.5 * (_elasticity_tensor[_qp](0, 0, 0, 0) -  _elasticity_tensor[_qp](0, 0, 1, 1));
    xm_temp = 9 *_xm_cali * k_b * _temperature[_qp] / (_burgers_vector_mag * _burgers_vector_mag * _burgers_vector_mag * sqrt(mu_f * eta_f));
	  
	if (_temperature[_qp] < _melting_temperature_high) { 
    
      if (xm_temp <= _xm_max) {
        if(xm_temp >= _xm) {
          xm = xm_temp; 
        } else {
          xm = _xm;
        }
      } else {
        xm = _xm_max;
      }
    } else {
      xm = _xm;
    }
  }

  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress = _tau[_qp][i] - _backstress[_qp][i];	  
	  
    if (MooseUtils::absoluteFuzzyEqual(effective_stress, 0.0)) {
		
      dslip_dtau[i] = 0.0;
      		
	} else {
		
	  stress_ratio = std::abs(effective_stress / _slip_resistance[_qp][i]);

      dslip_dtau[i] = _ao / xm *
                      std::pow(stress_ratio, 1.0 / xm - 1.0) /
                      _slip_resistance[_qp][i]
                    + creep_ao / _creep_xm *
                      std::pow(stress_ratio, 1.0 / _creep_xm - 1.0) /
                      _slip_resistance[_qp][i];
    }
  }
}

bool
CrystalPlasticityDislocationUpdateAluminum::areConstitutiveStateVariablesConverged()
{
  return isConstitutiveStateVariableConverged(_rho_ssd[_qp],
                                              _rho_ssd_before_update,
                                              _previous_substep_rho_ssd,
                                              _rho_tol);
}

void
CrystalPlasticityDislocationUpdateAluminum::updateSubstepConstitutiveVariableValues()
{
  // Update temporary variables at the end of the substep
  _previous_substep_rho_ssd = _rho_ssd[_qp];
  _previous_substep_rho_gnd_edge = _rho_gnd_edge[_qp];
  _previous_substep_rho_gnd_screw = _rho_gnd_screw[_qp];
  _previous_substep_backstress = _backstress[_qp];
}

void
CrystalPlasticityDislocationUpdateAluminum::cacheStateVariablesBeforeUpdate()
{
  _rho_ssd_before_update = _rho_ssd[_qp];
  _rho_gnd_edge_before_update = _rho_gnd_edge[_qp];
  _rho_gnd_screw_before_update = _rho_gnd_screw[_qp];
  _backstress_before_update = _backstress[_qp];
}

void
CrystalPlasticityDislocationUpdateAluminum::calculateStateVariableEvolutionRateComponent()
{
  Real rho_sum;

  // SSD dislocation density increment
  for (const auto i : make_range(_number_slip_systems))
  {
    
    rho_sum = _rho_ssd[_qp][i] + std::abs(_rho_gnd_edge[_qp][i]) + std::abs(_rho_gnd_screw[_qp][i]);

    // Multiplication and annihilation
	// note that _slip_increment here is the rate
	// and the rate equation gets multiplied by time step in updateStateVariables
    _rho_ssd_increment[i] = _k_0 * sqrt(rho_sum) - 2 * _y_c * _rho_ssd[_qp][i];
    _rho_ssd_increment[i] *= std::abs(_slip_increment[_qp][i]) / _burgers_vector_mag;

  }
  
  // GND dislocation density increment
  for (const auto i : make_range(_number_slip_systems)) 
  {

    _rho_gnd_edge_increment[i] = (-1.0) * _dslip_increment_dedge[_qp](i) / _burgers_vector_mag;
    _rho_gnd_screw_increment[i] = _dslip_increment_dscrew[_qp](i) / _burgers_vector_mag;
    
  }
  
  // backstress increment
  ArmstrongFrederickBackstressUpdate();
}

// Armstrong-Frederick update of the backstress
void
CrystalPlasticityDislocationUpdateAluminum::ArmstrongFrederickBackstressUpdate()
{
  for (const auto i : make_range(_number_slip_systems)) 
  {
    _backstress_increment[i] = _h * _slip_increment[_qp][i];
    _backstress_increment[i] -= _h_D * _backstress[_qp][i] * std::abs(_slip_increment[_qp][i]);  
  }
}

bool
CrystalPlasticityDislocationUpdateAluminum::updateStateVariables()
{
  // Now perform the check to see if the slip system should be updated
  // SSD
  for (const auto i : make_range(_number_slip_systems))
  { 
    _rho_ssd_increment[i] *= _substep_dt;

    // force positive SSD density
    if (_previous_substep_rho_ssd[i] < _zero_tol && _rho_ssd_increment[i] < 0.0)
      _rho_ssd[_qp][i] = _previous_substep_rho_ssd[i];
    else
      _rho_ssd[_qp][i] = _previous_substep_rho_ssd[i] + _rho_ssd_increment[i];

    if (_rho_ssd[_qp][i] < 0.0)
      return false;
  }
  
  // GND edge: can be both positive or negative
  for (const auto i : make_range(_number_slip_systems))
  { 
    _rho_gnd_edge_increment[i] *= _substep_dt;
    _rho_gnd_edge[_qp][i] = _previous_substep_rho_gnd_edge[i] + _rho_gnd_edge_increment[i];
  }
  
  // GND screw: can be both positive or negative
  for (const auto i : make_range(_number_slip_systems))
  { 
    _rho_gnd_screw_increment[i] *= _substep_dt;
    _rho_gnd_screw[_qp][i] = _previous_substep_rho_gnd_screw[i] + _rho_gnd_screw_increment[i];
  }
  
  // Backstress: can be both positive or negative
  for (const auto i : make_range(_number_slip_systems))
  { 
    _backstress_increment[i] *= _substep_dt;
    _backstress[_qp][i] = _previous_substep_backstress[i] + _backstress_increment[i];
  }
  
  return true;
}
