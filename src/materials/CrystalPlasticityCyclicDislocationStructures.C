// Marco Pelegatti
// Università di Udine
// Nicolò Grilli
// Università di Bristol
// 23 Marzo 2024

#include "CrystalPlasticityCyclicDislocationStructures.h"
#include "libmesh/int_range.h"
#include "libmesh/utility.h"
#include <cmath>
#include "Function.h"
#include "RotationTensor.h"

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
  params.addParam<Real>("shear_modulus", 86000.0, "Shear modulus in Taylor hardening law G");
  params.addParam<Real>("s_0", 256.0, "Thermal slip resistance");
  params.addParam<Real>("nu", 0.3, "Poisson's ratio for backstress calculation by Eshelby's inclusion");
  params.addParam<Real>("K_struct",3.0,"Constant of similitude for dislocation substructure");
  params.addParam<Real>("tau_0", 80.0, "Resolved shear stress at initial yield");
  params.addParam<Real>("k_c",2.0,"Coefficient K in channel dislocations evolution, representing accumulation rate");
  params.addParam<Real>("y_s",0.013,"Critical annihilation diameter for screw dislocations");
  params.addParam<Real>("f_0",0.5,"Initial dislocation walls volume fraction");
  params.addParam<Real>("f_inf",0.25,"Saturated dislocation walls volume fraction");
  params.addParam<Real>("k_f",2.0,"Rate constant for the evolution of the dislocation walls volume fraction");
  params.addParam<Real>("eta_0",50,"Initial Max/Min axis length ratio of the dislocation substructure");
  params.addParam<Real>("eta_inf",1,"Saturated Max/Min axis length ratio of the dislocation substructure");
  params.addParam<Real>("X_norm",400,"Normalization constant for the evolution of the dislocation substructure");
  params.addParam<Real>("init_d_struct",10.0,"Initial characteristic dislocation substructure length");
  params.addParam<Real>("init_rho_c",0.01,"Initial channel dislocation density");
  params.addParam<Real>("init_rho_w",0.01,"Initial wall dislocation density");
  params.addParam<Real>("init_rho_PSB",0.01,"Initial PSB dislocation density");
  params.addParam<Real>("k_w",2.0,"Coefficient K in wall dislocations evolution, representing accumulation rate");
  params.addParam<Real>("y_e",0.003,"Critical annihilation diameter for edge dislocations");
  params.addRequiredParam<std::vector<Real>>("B_ii", "Initial macroscopic backstress tensor components");
  params.addParam<Real>("f_PSB_0",0.0,"Initial PSB fraction");
  params.addParam<Real>("f_PSB_inf",0.2,"PSB fraction at stabilization");
  params.addParam<Real>("k_PSB",0.11,"Increasing rate of PSB fraction");
  params.addParam<Real>("epsilon_p_eff_cum_PSB",0.06,"Critical accumulated plastic strain to develop PSBs");
  params.addParam<Real>("eta_PSB",20,"Max/Min axis length ratio of the PSB");
  params.addParam<Real>("f_w_PSB",0.42,"PSB dislocation walls volume fraction");
  params.addParam<Real>("init_d_struct_PSB",1.0,"Characteristic PSB length");
  params.addParam<Real>("k_c_PSB",1.0,"Coefficient K in PSB dislocations evolution, representing accumulation rate");
  params.addParam<Real>("y_PSB",0.015,"Critical annihilation diameter for dislocations in PSBs");
  
  params.addParam<UserObjectName>("read_init_d",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element value "
                                  "of dislocation structure size");
  
  params.addParam<UserObjectName>("read_prop_user_object",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element "
                                  "specific property values from file");
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
	_s_0(getParam<Real>("s_0")),
	_nu(getParam<Real>("nu")),
	_K_struct(getParam<Real>("K_struct")),
	_tau_0(getParam<Real>("tau_0")),
	_k_c(getParam<Real>("k_c")),
	_y_s(getParam<Real>("y_s")),
	_f_0(getParam<Real>("f_0")),
	_f_inf(getParam<Real>("f_inf")),
	_k_f(getParam<Real>("k_f")),
	_eta_0(getParam<Real>("eta_0")),
	_eta_inf(getParam<Real>("eta_inf")),
	_X_norm(getParam<Real>("X_norm")),
	_init_d_struct(getParam<Real>("init_d_struct")),
	
	// Initial values of the state variables
    _init_rho_c(getParam<Real>("init_rho_c")),
	_init_rho_w(getParam<Real>("init_rho_w")),
	_init_rho_PSB(getParam<Real>("init_rho_PSB")),
	
	_k_w(getParam<Real>("k_w")),
	_y_e(getParam<Real>("y_e")),
	
	// PSBs variables and parameters
	_f_PSB_0(getParam<Real>("f_PSB_0")),
	_f_PSB_inf(getParam<Real>("f_PSB_inf")),
	_k_PSB(getParam<Real>("k_PSB")),
	_epsilon_p_eff_cum_PSB(getParam<Real>("epsilon_p_eff_cum_PSB")),
	_eta_PSB(getParam<Real>("eta_PSB")),
	_f_w_PSB(getParam<Real>("f_w_PSB")),
	_init_d_struct_PSB(getParam<Real>("init_d_struct_PSB")),
	_k_c_PSB(getParam<Real>("k_c_PSB")),
	_y_PSB(getParam<Real>("y_PSB")),
	
	// State variables of the dislocation model
    _rho_c(declareProperty<std::vector<Real>>("rho_c")),
    _rho_c_old(getMaterialPropertyOld<std::vector<Real>>("rho_c")),
	_rho_w(declareProperty<std::vector<Real>>("rho_w")),
    _rho_w_old(getMaterialPropertyOld<std::vector<Real>>("rho_w")),
	_rho_PSB(declareProperty<std::vector<Real>>("rho_PSB")),
    _rho_PSB_old(getMaterialPropertyOld<std::vector<Real>>("rho_PSB")),
    
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
    _d_struct_old(getMaterialPropertyOld<Real>("d_struct")),
	
	// Mean glide distance for dislocations in the channel phase
    _l_c(declareProperty<Real>("l_c")),    
    
	// Characteristic dislocation substructure length in PSB
    _d_struct_PSB(declareProperty<Real>("d_struct_PSB")),
	_d_struct_PSB_old(getMaterialPropertyOld<Real>("d_struct_PSB")),
	
	// Mean glide distance for dislocations in the PSB
    _l_PSB(declareProperty<Real>("l_PSB")),
	
	// PSB fraction
    _f_PSB(declareProperty<Real>("f_PSB")),
	_f_PSB_old(getMaterialPropertyOld<Real>("f_PSB")),
	
    // Backstress variables
    _backstress_c(declareProperty<std::vector<Real>>("backstress_c")),
    _backstress_c_old(getMaterialPropertyOld<std::vector<Real>>("backstress_c")),
	_backstress_w(declareProperty<std::vector<Real>>("backstress_w")),
    _backstress_w_old(getMaterialPropertyOld<std::vector<Real>>("backstress_w")),
	_backstress_PSB(declareProperty<std::vector<Real>>("backstress_PSB")),
    _backstress_PSB_old(getMaterialPropertyOld<std::vector<Real>>("backstress_PSB")),
	
	// Slip resistance in channel, wall and PSB
	_slip_resistance_c(declareProperty<std::vector<Real>>(_base_name + "slip_resistance_c")),
    _slip_resistance_c_old(getMaterialPropertyOld<std::vector<Real>>(_base_name + "slip_resistance_c")),
	_slip_resistance_w(declareProperty<std::vector<Real>>(_base_name + "slip_resistance_w")),
    _slip_resistance_w_old(getMaterialPropertyOld<std::vector<Real>>(_base_name + "slip_resistance_w")),
	_slip_resistance_PSB(declareProperty<std::vector<Real>>(_base_name + "slip_resistance_PSB")),
    _slip_resistance_PSB_old(getMaterialPropertyOld<std::vector<Real>>(_base_name + "slip_resistance_PSB")),
	
	// Slip increment inside the channel, the wall and PSB
	_slip_increment_c(declareProperty<std::vector<Real>>("slip_increment_c")),
	_slip_increment_w(declareProperty<std::vector<Real>>("slip_increment_w")),
	_slip_increment_PSB(declareProperty<std::vector<Real>>("slip_increment_PSB")),

    // increments of state variables
    _rho_c_increment(_number_slip_systems, 0.0),
	_rho_w_increment(_number_slip_systems, 0.0),
	_rho_PSB_increment(_number_slip_systems, 0.0),
    _backstress_c_increment(_number_slip_systems, 0.0),
	_backstress_w_increment(_number_slip_systems, 0.0),
	_backstress_PSB_increment(_number_slip_systems, 0.0),
	
	// resize local caching vectors used for substepping
	_previous_substep_rho_c(_number_slip_systems, 0.0),
	_previous_substep_rho_w(_number_slip_systems, 0.0),
	_previous_substep_rho_PSB(_number_slip_systems, 0.0),
	_previous_substep_backstress_c(_number_slip_systems, 0.0),
	_previous_substep_backstress_w(_number_slip_systems, 0.0),
	_previous_substep_backstress_PSB(_number_slip_systems, 0.0),
    _rho_c_before_update(_number_slip_systems, 0.0),
	_rho_w_before_update(_number_slip_systems, 0.0),
	_rho_PSB_before_update(_number_slip_systems, 0.0),
    _backstress_c_before_update(_number_slip_systems, 0.0),
	_backstress_w_before_update(_number_slip_systems, 0.0),
	_backstress_PSB_before_update(_number_slip_systems, 0.0),
	
	// Interaction matrix between slip systems
	_A_int(_number_slip_systems, _number_slip_systems),
	
	// Initial macroscopic backstress tensor components
	_B_ii(getParam<std::vector<Real>>("B_ii")),
	_B_0(_B_ii[0], _B_ii[1], _B_ii[2], _B_ii[3], _B_ii[4], _B_ii[5]),
	
	// Element property read user object used to read initial dislocation structure size
    _read_init_d(isParamValid("read_init_d")
                               ? &getUserObject<PropertyReadFile>("read_init_d")
                               : nullptr),
	
	// Element property read user object used to read in Euler angles
    _read_prop_user_object(isParamValid("read_prop_user_object")
                               ? &getUserObject<PropertyReadFile>("read_prop_user_object")
                               : nullptr)
{
}

void
CrystalPlasticityCyclicDislocationStructures::initQpStatefulProperties()
{
  // Slip resistance is resized here (not anymore)
  CrystalPlasticityDislocationUpdateBase::initQpStatefulProperties();
  
  Real taylor_hardening;

  // Initialize the dislocation density size
  _rho_c[_qp].resize(_number_slip_systems);
  _rho_w[_qp].resize(_number_slip_systems);
  _rho_PSB[_qp].resize(_number_slip_systems);
  
  // Initialize the backstress size
  _backstress_c[_qp].resize(_number_slip_systems);
  _backstress_w[_qp].resize(_number_slip_systems);
  _backstress_PSB[_qp].resize(_number_slip_systems);
  
  // Initialize slip resistance inside channel, wall and PSB
  _slip_resistance_c[_qp].resize(_number_slip_systems);
  _slip_resistance_w[_qp].resize(_number_slip_systems);
  _slip_resistance_PSB[_qp].resize(_number_slip_systems);
  
  // Initialize slip increment inside channel, wall and PSB
  _slip_increment_c[_qp].resize(_number_slip_systems);
  _slip_increment_w[_qp].resize(_number_slip_systems);
  _slip_increment_PSB[_qp].resize(_number_slip_systems);
  
  // Initialize instantaneous plastic deformation tangent at the slip system level
  _dslip_dtau[_qp].resize(_number_slip_systems);
  
  // Initialize walls fraction
  _f_w[_qp] = _f_0;
  
  // Initialize max/min axis length ratio of the dislocation substructure
  _eta[_qp] = _eta_0;
  
  // Initialize characteristic dislocation substructure length
  
  if (_read_init_d) { // Read initial characteristic dislocation substructure length from file
  
    _d_struct[_qp] = _read_init_d->getData(_current_elem, 0);
    
    } else { // Initialize uniform characteristic dislocation substructure length
    
    _d_struct[_qp] = _init_d_struct;
	// _d_struct_PSB[_qp] = _init_d_struct_PSB;
    
  }
  
  // Initialize characteristic dislocation substructure length in PSB
  
  _d_struct_PSB[_qp] = _d_struct[_qp];
  
  // Initialize mean glide distance for dislocations in the channel phase
  _l_c[_qp] = _eta[_qp] * _d_struct[_qp];
  
  // Initialize mean glide distance for dislocations in the PSB
  _l_PSB[_qp] = _eta_PSB * _d_struct_PSB[_qp];
  
  // Initialize PSB fraction
  _f_PSB[_qp] = _f_PSB_0;
  
  assignEulerAngles();
  
  calculateFlowDirection(_crysrot);
  
  // Initialize dislocation densities and backstress
  for (const auto i : make_range(_number_slip_systems))
  {
    _rho_c[_qp][i] = _init_rho_c;
	_rho_w[_qp][i] = _init_rho_w;
	_rho_PSB[_qp][i] = _init_rho_PSB;
	
	_backstress_c[_qp][i] = _f_w[_qp] * _B_0.doubleContraction(_flow_direction[_qp][i]);
    _backstress_w[_qp][i] = - ( 1.0 - _f_w[_qp] ) * _B_0.doubleContraction(_flow_direction[_qp][i]);
	_backstress_PSB[_qp][i] = _f_PSB[_qp] * _B_0.doubleContraction(_flow_direction[_qp][i]);
  }
  
  initializeInteractionMatrix();

  // Initialize value of the slip resistance
  // as a function of the dislocation density
  for (const auto i : make_range(_number_slip_systems))
  {
    // Add Peierls stress
    _slip_resistance_c[_qp][i] = _s_0;
	
	_slip_resistance_c[_qp][i] += (_shear_modulus * _burgers_vector_mag
	                           * std::sqrt(_A_int(i,i) * _rho_c[_qp][i]));
							   
	// Add Peierls stress
	_slip_resistance_w[_qp][i] = _s_0;
	
	taylor_hardening = 0.0;
	
	for (const auto j : make_range(_number_slip_systems))
    {
      taylor_hardening += _A_int(i,j) * _rho_w[_qp][j];
	}
	
    _slip_resistance_w[_qp][i] += _shear_modulus * _burgers_vector_mag
	                           * std::sqrt(taylor_hardening);
	
	_slip_resistance_PSB[_qp][i] = _s_0;
	
	_slip_resistance_PSB[_qp][i] += (_shear_modulus * _burgers_vector_mag
	                           * std::sqrt(_A_int(i,i) * _rho_PSB[_qp][i]));
  }

  // initialize slip increment
  for (const auto i : make_range(_number_slip_systems))
  {
	_slip_increment_c[_qp][i] = 0.0;
	_slip_increment_w[_qp][i] = 0.0;
	_slip_increment_PSB[_qp][i] = 0.0;
	
	_slip_increment[_qp][i] = 0.0;
  }
}

void
CrystalPlasticityCyclicDislocationStructures::initializeInteractionMatrix()
{
  Real G0 = 0.122;
  Real G1 = 0.122;
  Real G2 = 0.625;
  Real G3 = 0.137;
  Real G4 = 0.07;
  Real G5 = 0.122;
  
  Real A_int[_number_slip_systems][_number_slip_systems] = {
  {G0,G1,G1,G2,G3,G3,G4,G3,G5,G4,G5,G3},
  {G1,G0,G1,G3,G4,G5,G3,G2,G3,G5,G4,G3},
  {G1,G1,G0,G3,G5,G4,G5,G3,G4,G3,G3,G2},
  {G2,G3,G3,G0,G1,G1,G4,G5,G3,G4,G3,G5},
  {G3,G4,G5,G1,G0,G1,G5,G4,G3,G3,G2,G3},
  {G3,G5,G4,G1,G1,G0,G3,G3,G2,G5,G3,G4},
  {G4,G3,G5,G4,G5,G3,G0,G1,G1,G2,G3,G3},
  {G3,G2,G3,G5,G4,G3,G1,G0,G1,G3,G4,G5},
  {G5,G3,G4,G3,G3,G2,G1,G1,G0,G3,G5,G4},
  {G4,G5,G3,G4,G3,G5,G2,G3,G3,G0,G1,G1},
  {G5,G4,G3,G3,G2,G3,G3,G4,G5,G1,G0,G1},
  {G3,G3,G2,G5,G3,G4,G3,G5,G4,G1,G1,G0}};
  
  for (const auto i : make_range(_number_slip_systems))
    for (const auto j : make_range(_number_slip_systems))
      _A_int(i,j) = A_int[i][j];
}

void
CrystalPlasticityCyclicDislocationStructures::assignEulerAngles()
{ 
  for (const auto i : make_range(LIBMESH_DIM)) {
    _Euler_angles(i) = _read_prop_user_object->getData(_current_elem, i);	  
  }

  RotationTensor R(_Euler_angles);
  _crysrot = R.transpose();
}

void
CrystalPlasticityCyclicDislocationStructures::setInitialConstitutiveVariableValues()
{
  // Initialize state variables with the value at the previous time step
  _rho_c[_qp] = _rho_c_old[_qp];
  _previous_substep_rho_c = _rho_c_old[_qp];
  _rho_w[_qp] = _rho_w_old[_qp];
  _previous_substep_rho_w = _rho_w_old[_qp];
  _rho_PSB[_qp] = _rho_PSB_old[_qp];
  _previous_substep_rho_PSB = _rho_PSB_old[_qp];
  _backstress_c[_qp] = _backstress_c_old[_qp];
  _previous_substep_backstress_c = _backstress_c_old[_qp];
  _backstress_w[_qp] = _backstress_w_old[_qp];
  _previous_substep_backstress_w = _backstress_w_old[_qp];
  _backstress_PSB[_qp] = _backstress_PSB_old[_qp];
  _previous_substep_backstress_PSB = _backstress_PSB_old[_qp];
}

void
CrystalPlasticityCyclicDislocationStructures::setSubstepConstitutiveVariableValues()
{
  // Inialize state variable of the next substep
  // with the value at the previous substep
  _rho_c[_qp] = _previous_substep_rho_c;
  _rho_w[_qp] = _previous_substep_rho_w;
  _rho_PSB[_qp] = _previous_substep_rho_PSB;
  _backstress_c[_qp] = _previous_substep_backstress_c;
  _backstress_w[_qp] = _previous_substep_backstress_w;
  _backstress_PSB[_qp] = _previous_substep_backstress_PSB;
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
  Real stress_ratio_c;
  Real stress_ratio_w;
  Real stress_ratio_PSB;
  
  // Difference between RSS and backstress
  // temporary variable for each slip system
  Real effective_stress_c;
  Real effective_stress_w;
  Real effective_stress_PSB;
  
  // Slip increment within the matrix
  Real slip_increment_matrix;
  
  // Slip rate inside channel
  
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress_c = _tau[_qp][i] - _backstress_c[_qp][i];
    
    stress_ratio_c = std::abs(effective_stress_c / _slip_resistance_c[_qp][i]);
    
    _slip_increment_c[_qp][i] =
        _ao * std::pow(stress_ratio_c, 1.0 / _xm);
      
    if (effective_stress_c < 0.0)
      _slip_increment_c[_qp][i] *= -1.0;

  }
  
  // Slip rate inside wall
  
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress_w = _tau[_qp][i] - _backstress_w[_qp][i];
    
    stress_ratio_w = std::abs(effective_stress_w / _slip_resistance_w[_qp][i]);
    
    _slip_increment_w[_qp][i] =
        _ao * std::pow(stress_ratio_w, 1.0 / _xm);
      
    if (effective_stress_w < 0.0)
      _slip_increment_w[_qp][i] *= -1.0;

  }
  
  // Slip rate inside PSB
  
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress_PSB = _tau[_qp][i] - _backstress_PSB[_qp][i];
    
    stress_ratio_PSB = std::abs(effective_stress_PSB / _slip_resistance_PSB[_qp][i]);
    
    _slip_increment_PSB[_qp][i] =
        _ao * std::pow(stress_ratio_PSB, 1.0 / _xm);
      
    if (effective_stress_PSB < 0.0)
      _slip_increment_PSB[_qp][i] *= -1.0;

  }
  
  // Total slip increment
  
  for (const auto i : make_range(_number_slip_systems))
  {
    slip_increment_matrix = (1.0 - _f_w[_qp]) * _slip_increment_c[_qp][i] + _f_w[_qp] * _slip_increment_w[_qp][i];
	
	_slip_increment[_qp][i] = (1.0 - _f_PSB[_qp]) * slip_increment_matrix + _f_PSB[_qp] * _slip_increment_PSB[_qp][i];

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
    _slip_resistance_c[_qp][i] = _s_0;
	
	_slip_resistance_c[_qp][i] += (_shear_modulus * _burgers_vector_mag
	                           * std::sqrt(_A_int(i,i) * _rho_c[_qp][i]));
							   
	// Add Peierls stress
	_slip_resistance_w[_qp][i] = _s_0;
	
	taylor_hardening = 0.0;
	
	for (const auto j : make_range(_number_slip_systems))
    {
      taylor_hardening += _A_int(i,j) * _rho_w[_qp][j];
	}
	
    _slip_resistance_w[_qp][i] += _shear_modulus * _burgers_vector_mag
	                           * std::sqrt(taylor_hardening);
	
	_slip_resistance_PSB[_qp][i] = _s_0;
	
	_slip_resistance_PSB[_qp][i] += (_shear_modulus * _burgers_vector_mag
	                           * std::sqrt(_A_int(i,i) * _rho_PSB[_qp][i]));
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
  Real stress_ratio_c;
  Real stress_ratio_w;
  Real stress_ratio_PSB;
  
  // Difference between RSS and backstress
  // temporary variable for each slip system
  Real effective_stress_c;
  Real effective_stress_w;
  Real effective_stress_PSB;
  
  // Derivative slip rates
  // temporary variable for each slip system
  std::vector<Real> dslip_dtau_c(_number_slip_systems, 0.0);
  std::vector<Real> dslip_dtau_w(_number_slip_systems, 0.0);
  std::vector<Real> dslip_dtau_matrix(_number_slip_systems, 0.0);
  std::vector<Real> dslip_dtau_PSB(_number_slip_systems, 0.0);
  
  // Derivative slip rate inside channel
	
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress_c = _tau[_qp][i] - _backstress_c[_qp][i];	  
	  
    if (MooseUtils::absoluteFuzzyEqual(effective_stress_c, 0.0)) {
		
      dslip_dtau_c[i] = 0.0;
      		
	} else {
		
	  stress_ratio_c = std::abs(effective_stress_c / _slip_resistance_c[_qp][i]);

      dslip_dtau_c[i] = _ao / _xm *
                      std::pow(stress_ratio_c, 1.0 / _xm - 1.0) /
                      _slip_resistance_c[_qp][i];		
	}
  }
  
  // Derivative slip rate inside wall
	
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress_w = _tau[_qp][i] - _backstress_w[_qp][i];	  
	  
    if (MooseUtils::absoluteFuzzyEqual(effective_stress_w, 0.0)) {
		
      dslip_dtau_w[i] = 0.0;
      		
	} else {
		
	  stress_ratio_w = std::abs(effective_stress_w / _slip_resistance_w[_qp][i]);

      dslip_dtau_w[i] = _ao / _xm *
                      std::pow(stress_ratio_w, 1.0 / _xm - 1.0) /
                      _slip_resistance_w[_qp][i];		
	}
  }
  
  // Derivative slip rate inside PSB
	
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress_PSB = _tau[_qp][i] - _backstress_PSB[_qp][i];	  
	  
    if (MooseUtils::absoluteFuzzyEqual(effective_stress_PSB, 0.0)) {
		
      dslip_dtau_PSB[i] = 0.0;
      		
	} else {
		
	  stress_ratio_PSB = std::abs(effective_stress_PSB / _slip_resistance_PSB[_qp][i]);

      dslip_dtau_PSB[i] = _ao / _xm *
                      std::pow(stress_ratio_PSB, 1.0 / _xm - 1.0) /
                      _slip_resistance_PSB[_qp][i];		
	}
  }
  
  // Derivative total slip rate
	
  for (const auto i : make_range(_number_slip_systems))
  {
	dslip_dtau_matrix[i] = (1.0 - _f_w[_qp]) * dslip_dtau_c[i] + _f_w[_qp] * dslip_dtau_w[i];
	
	dslip_dtau[i] = (1.0 - _f_PSB[_qp]) * dslip_dtau_matrix[i] + _f_PSB[_qp] * dslip_dtau_PSB[i];
	
	_dslip_dtau[_qp][i] = dslip_dtau[i];
  }
}

void
CrystalPlasticityCyclicDislocationStructures::updateSubstepConstitutiveVariableValues()
{
  // Update temporary variables at the end of the substep
  _previous_substep_rho_c = _rho_c[_qp];
  _previous_substep_rho_w = _rho_w[_qp];
  _previous_substep_rho_PSB = _rho_PSB[_qp];
  _previous_substep_backstress_c = _backstress_c[_qp];
  _previous_substep_backstress_w = _backstress_w[_qp];
  _previous_substep_backstress_PSB = _backstress_PSB[_qp];
}

void
CrystalPlasticityCyclicDislocationStructures::cacheStateVariablesBeforeUpdate()
{
  _rho_c_before_update = _rho_c[_qp];
  _rho_w_before_update = _rho_w[_qp];
  _rho_PSB_before_update = _rho_PSB[_qp];
  _backstress_c_before_update = _backstress_c[_qp];
  _backstress_w_before_update = _backstress_w[_qp];
  _backstress_PSB_before_update = _backstress_PSB[_qp];
}

// Update channel dislocation density based on (27) - (30)
void
CrystalPlasticityCyclicDislocationStructures::calculateStateVariableEvolutionRateComponent()
{
  // Update substructure variables
  calculateWallVolumeFraction();
  calculateSubstructureParameter();
  calculateSubstructureSize();
  calculatePSBFraction();
  
  // Slip increment within the matrix
  Real slip_increment_matrix;

  // channel and wall dislocation density increment
  for (const auto i : make_range(_number_slip_systems))
  {
    // Multiplication and annihilation
	// note that _slip_increment here is the rate
	// and the rate equation gets multiplied by time step in updateStateVariables
    _rho_c_increment[i] = _k_c / _l_c[_qp] - 2.0 * _y_s * _rho_c[_qp][i];
    _rho_c_increment[i] *= std::abs(_slip_increment_c[_qp][i]) / _burgers_vector_mag;
	
	slip_increment_matrix = (1.0 - _f_w[_qp]) * _slip_increment_c[_qp][i] + _f_w[_qp] * _slip_increment_w[_qp][i];
	
	_rho_w_increment[i] = ( _k_w * std::abs(_slip_increment_c[_qp][i]) ) / ( _burgers_vector_mag * _l_c[_qp] );
    _rho_w_increment[i] -= ( 2.0 * _y_e * _rho_w[_qp][i] * std::abs(slip_increment_matrix) ) / _burgers_vector_mag;
	
	if (_epsilon_p_eff_cum[_qp] > _epsilon_p_eff_cum_PSB) {
  
      _rho_PSB_increment[i] = _k_c_PSB / _l_PSB[_qp] - 2.0 * _y_PSB * _rho_PSB[_qp][i];
	  _rho_PSB_increment[i] *= std::abs(_slip_increment_PSB[_qp][i]) / _burgers_vector_mag;
  
    } else {
	  
      _rho_PSB_increment[i] = _rho_c_increment[i];
	  
    }
  }
  
  // backstress increment
  BackstressUpdate();
}

// Backstress update based on (15) e (16)
void
CrystalPlasticityCyclicDislocationStructures::BackstressUpdate()
{
  // slip system independent shape factor
  Real K_shape;
  Real K_shape_PSB;

  // degree of plastic accommodation at the interface of the two phases
  Real f_accom;
  
  // modified Poisson's ratio
  Real nu_p;

  // S1212 is a component of the Eshelby tensor
  Real S_1212;
  Real S_1212_PSB;
  
  // modified shear modulus
  Real mu_m;
	
  // Function that depends on the shape of the channel phase
  Real f_eta;
  Real f_eta_PSB;
  
  K_shape = _eta[_qp] * std::sqrt( std::pow(_eta[_qp], 2.0) - 1.0) - std::acosh(_eta[_qp]);
  K_shape *= 2.0 * libMesh::pi * _eta[_qp];
  K_shape /= std::sqrt(std::pow(std::pow(_eta[_qp], 2.0) - 1.0, 3.0));
  
  K_shape_PSB = _eta_PSB * std::sqrt( std::pow(_eta_PSB, 2.0) - 1.0) - std::acosh(_eta_PSB);
  K_shape_PSB *= 2.0 * libMesh::pi * _eta_PSB;
  K_shape_PSB /= std::sqrt(std::pow(std::pow(_eta_PSB, 2.0) - 1.0, 3.0));
  
  for (const auto i : make_range(_number_slip_systems)) 
  {
    f_accom = _dslip_dtau[_qp][i] * _substep_dt / 2.0;
    
    nu_p = _nu + (2.0/3.0) * (1.0 + _nu) * _shear_modulus * f_accom;
	nu_p /= 1.0 + (4.0/3.0) * (1.0 + _nu) * _shear_modulus * f_accom;
	
	S_1212 = std::pow(_eta[_qp], 2.0) - 1.75 - 2.0 * nu_p * std::pow(_eta[_qp], 2.0) + nu_p;
	S_1212 *= K_shape;
	S_1212 += libMesh::pi * std::pow(_eta[_qp], 2.0);
	S_1212 /= 8.0 * libMesh::pi * (1.0 - nu_p) * (std::pow(_eta[_qp], 2.0) - 1.0);
	
	mu_m = _shear_modulus / (1.0 + 2.0 * _shear_modulus * f_accom);
	
	f_eta = mu_m * (1.0 - 2.0 * S_1212) / S_1212;
	f_eta /= 1.0 + (mu_m / _shear_modulus) * ((1.0 - 2.0 * S_1212) / (2.0 * S_1212));

    _backstress_c_increment[i] = _f_w[_qp] * f_eta * ( _slip_increment_c[_qp][i] - _slip_increment_w[_qp][i] );
	_backstress_w_increment[i] = - ( 1.0 - _f_w[_qp] ) * f_eta * ( _slip_increment_c[_qp][i] - _slip_increment_w[_qp][i] );
	
	S_1212_PSB = std::pow(_eta_PSB, 2.0) - 1.75 - 2.0 * nu_p * std::pow(_eta_PSB, 2.0) + nu_p;
	S_1212_PSB *= K_shape_PSB;
	S_1212_PSB += libMesh::pi * std::pow(_eta_PSB, 2.0);
	S_1212_PSB /= 8.0 * libMesh::pi * (1.0 - nu_p) * (std::pow(_eta_PSB, 2.0) - 1.0);
	
	f_eta_PSB = mu_m * (1.0 - 2.0 * S_1212_PSB) / S_1212_PSB;
	f_eta_PSB /= 1.0 + (mu_m / _shear_modulus) * ((1.0 - 2.0 * S_1212_PSB) / (2.0 * S_1212_PSB));
	
	_backstress_PSB_increment[i] = ( _f_w_PSB / ( 1.0 - _f_w_PSB ) ) * f_eta_PSB * _slip_increment_PSB[_qp][i];
  }
}

bool
CrystalPlasticityCyclicDislocationStructures::updateStateVariables()
{		
  for (const auto i : make_range(_number_slip_systems))
  { 
    _rho_c_increment[i] *= _substep_dt;
	_rho_w_increment[i] *= _substep_dt;
	_rho_PSB_increment[i] *= _substep_dt;

    // force positive channel dislocation density
    if (_previous_substep_rho_c[i] < _zero_tol && _rho_c_increment[i] < 0.0)
      _rho_c[_qp][i] = _previous_substep_rho_c[i];
    else
      _rho_c[_qp][i] = _previous_substep_rho_c[i] + _rho_c_increment[i];

    if (_rho_c[_qp][i] < 0.0)
      return false;
  
    // force positive wall dislocation density
    if (_previous_substep_rho_w[i] < _zero_tol && _rho_w_increment[i] < 0.0)
      _rho_w[_qp][i] = _previous_substep_rho_w[i];
    else
      _rho_w[_qp][i] = _previous_substep_rho_w[i] + _rho_w_increment[i];

    if (_rho_w[_qp][i] < 0.0)
      return false;
  
    // force positive PSB dislocation density
    if (_previous_substep_rho_PSB[i] < _zero_tol && _rho_PSB_increment[i] < 0.0)
      _rho_PSB[_qp][i] = _previous_substep_rho_PSB[i];
    else
      _rho_PSB[_qp][i] = _previous_substep_rho_PSB[i] + _rho_PSB_increment[i];

    if (_rho_PSB[_qp][i] < 0.0)
      return false;
  }
  
  // Backstress: can be both positive or negative
  for (const auto i : make_range(_number_slip_systems))
  { 
    _backstress_c_increment[i] *= _substep_dt;
    _backstress_c[_qp][i] = _previous_substep_backstress_c[i] + _backstress_c_increment[i];
	_backstress_w_increment[i] *= _substep_dt;
    _backstress_w[_qp][i] = _previous_substep_backstress_w[i] + _backstress_w_increment[i];
	_backstress_PSB_increment[i] *= _substep_dt;
    _backstress_PSB[_qp][i] = _previous_substep_backstress_PSB[i] + _backstress_PSB_increment[i];
  }
  
  return true;
}

// Calculate wall volume fraction based on (44)
void
CrystalPlasticityCyclicDislocationStructures::calculateWallVolumeFraction()
{
  _f_w[_qp] = _f_inf + (_f_0 - _f_inf) * std::exp(-_epsilon_p_eff_cum[_qp] / _k_f);
}

// Calculate substructure parameter eta based on (41), (42), (43)
void
CrystalPlasticityCyclicDislocationStructures::calculateSubstructureParameter()
{
  // Sum of channel dislocation density on all slip systems
  Real sum_rho_c = 0.0;
  
  // Maximum channel dislocation density among slip systems
  Real max_rho_c = 0.0;
  
  // k_p parameter for evolution rate of eta
  Real k_p;
	
  for (const auto i : make_range(_number_slip_systems))
  {
    sum_rho_c += _rho_c[_qp][i];
    
    if (_rho_c[_qp][i] > max_rho_c) {
		
      max_rho_c = _rho_c[_qp][i];
      
	}
  }
  
  k_p = max_rho_c / sum_rho_c;
  
  _eta[_qp] = _eta_inf + (_eta_0 - _eta_inf) * std::exp(-_X_norm * _epsilon_p_eff_cum[_qp] / k_p);
	
}

// Calculate substructure size based on (40)
void
CrystalPlasticityCyclicDislocationStructures::calculateSubstructureSize()
{
  // Max RSS among the slip systems	
  Real max_abs_tau = 0.0; 
  
  // Prefactor of _d_struct
  Real d_struct_prefactor = _K_struct * _shear_modulus * _burgers_vector_mag;
  
  for (const auto i : make_range(_number_slip_systems))
  {
    if (std::abs(_tau[_qp][i]) > max_abs_tau) {
		
		max_abs_tau = std::abs(_tau[_qp][i]);
		
	}
  }
  
  if (max_abs_tau > d_struct_prefactor / _d_struct_old[_qp] + _tau_0) { // _d_struct can only decrease
  
    _d_struct[_qp] = d_struct_prefactor / (max_abs_tau - _tau_0);
  
  } else {
	  
    _d_struct[_qp] = _d_struct_old[_qp];
	  
  }

  _l_c[_qp] = _eta[_qp] * _d_struct[_qp];
  
  _l_PSB[_qp] = _eta_PSB * _d_struct_PSB[_qp];
}

// Calculate PSB fraction
void
CrystalPlasticityCyclicDislocationStructures::calculatePSBFraction()
{  
  if (_epsilon_p_eff_cum[_qp] > _epsilon_p_eff_cum_PSB) {
  
	// _f_PSB[_qp] = _f_PSB_inf + (_f_PSB_0 - _f_PSB_inf) * std::exp( - _k_PSB * (_epsilon_p_eff_cum[_qp] - _epsilon_p_eff_cum_PSB) );
	
	_f_PSB[_qp] = _f_PSB_0 + _k_PSB * std::sqrt( _epsilon_p_eff_cum[_qp] - _epsilon_p_eff_cum_PSB );
  
  } else {
	  
    _f_PSB[_qp] = _f_PSB_old[_qp];
	  
  }
}
