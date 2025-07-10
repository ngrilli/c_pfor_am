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
                             "Includes slip and backstress. ");
  params.addParam<Real>("gamma_o", 0.004, "reference slip rate");
  params.addParam<Real>("m_exp", 0.1, "strain rate sensitivity exponent");
  params.addParam<bool>("cap_slip_increment", false, "Cap the absolute value of the slip increment "
                                                     "in one time step to _slip_incr_tol");												 
  params.addParam<Real>("burgers_vector_mag", 0.000256, "Magnitude of the Burgers vector");
  params.addParam<Real>("shear_modulus", 86000.0, "Shear modulus in Taylor hardening law G");
  params.addParam<Real>("s_0", 256.0, "Initial slip resistance");
  params.addParam<Real>("s_0c", 256.0, "Initial slip resistance in channel");
  params.addParam<Real>("s_0w", 256.0, "Initial slip resistance in wall");
  params.addParam<Real>("nu", 0.3, "Poisson's ratio for backstress calculation by Eshelby's inclusion");
  params.addParam<Real>("k_c",2.0,"Coefficient of accumulation rate of dislocations in channel phase");
  params.addParam<Real>("y_c",0.013,"Critical annihilation distance for dislocations in channel phase");
  params.addParam<Real>("k_w",2.0,"Coefficient of accumulation rate of dislocations in wall phase");
  params.addParam<Real>("y_w",0.003,"Critical annihilation distance for dislocations in wall phase");
  params.addParam<Real>("f_0",0.5,"Initial wall volume fraction");
  params.addParam<Real>("f_inf",0.25,"Saturated wall volume fraction");
  params.addParam<Real>("k_f",2.0,"Rate constant for evolution of wall volume fraction");
  params.addParam<Real>("eta_0",50,"Initial channel aspect ratio of dislocation structure");
  params.addParam<Real>("eta_inf",1,"Saturated channel aspect ratio of dislocation structure");
  params.addParam<std::vector<Real>>("eta_0_alpha", {}, "Initial channel aspect ratio of cellular structure at slip system level");
  params.addParam<Real>("X_norm",400,"Normalization constant for evolution of dislocation structure");
  params.addParam<Real>("K_struct",3.0,"Constant of similitude for dislocation structure");
  params.addParam<Real>("tau_0", 80.0, "Resolved shear stress at initial yield for similitude scaling law");
  params.addParam<Real>("init_d_struct",10.0,"Initial characteristic dislocation structure length");
  params.addParam<Real>("init_rho_c",0.01,"Initial channel dislocation density");
  params.addParam<Real>("init_rho_w",0.01,"Initial wall dislocation density");
  params.addRequiredParam<std::vector<Real>>("B_ii", "Initial macroscopic backstress tensor components");
  params.addParam<Real>("f_PSB_0",0.0,"Initial PSB volume fraction");
  params.addParam<Real>("f_PSB_inf",0.2,"Saturated PSB volume fraction");
  params.addParam<Real>("k_PSB",0.11,"Rate constant for evolution of PSB volume fraction");
  params.addParam<Real>("epsilon_p_eff_cum_PSB",0.06,"Critical cumulative effective plastic strain to develop PSB");
  params.addParam<Real>("k_c_PSB",1.0,"Coefficient of accumulation rate of dislocations in PSB channel phase");
  params.addParam<Real>("y_PSB",0.015,"Critical annihilation distance for dislocations in PSB phase");
  params.addParam<Real>("f_w_PSB",0.42,"PSB wall volume fraction");
  params.addParam<Real>("eta_PSB",20,"Channel aspect ratio of PSB");
  params.addParam<Real>("init_d_struct_PSB",1.0,"Initial characteristic PSB length");
  params.addParam<Real>("init_rho_PSB",0.01,"Initial PSB dislocation density");
  
  params.addParam<Real>("A_self",0.122,"Self interaction coefficient");
  params.addParam<Real>("A_copl",0.122,"Coplanar interaction coefficient");
  params.addParam<Real>("A_CS",0.625,"Collinear interaction coefficient");
  params.addParam<Real>("A_GJ",0.137,"Glissile junction interaction coefficient");
  params.addParam<Real>("A_HL",0.07,"Hirth lock interaction coefficient");
  params.addParam<Real>("A_LC",0.122,"Lomer-Cottrel lock interaction coefficient");
  
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
    _gamma_o(getParam<Real>("gamma_o")),
    _m_exp(getParam<Real>("m_exp")),
    _cap_slip_increment(getParam<bool>("cap_slip_increment")),
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")),
	_shear_modulus(getParam<Real>("shear_modulus")),
	_s_0(getParam<Real>("s_0")),
	_s_0c(getParam<Real>("s_0c")),
	_s_0w(getParam<Real>("s_0w")),
	_nu(getParam<Real>("nu")),
	_k_c(getParam<Real>("k_c")),
	_y_c(getParam<Real>("y_c")),
	_k_w(getParam<Real>("k_w")),
	_y_w(getParam<Real>("y_w")),
	_f_0(getParam<Real>("f_0")),
	_f_inf(getParam<Real>("f_inf")),
	_k_f(getParam<Real>("k_f")),
	_eta_0(getParam<Real>("eta_0")),
	_eta_inf(getParam<Real>("eta_inf")),
	_eta_0_alpha(getParam<std::vector<Real>>("eta_0_alpha")),
	_X_norm(getParam<Real>("X_norm")),
	_K_struct(getParam<Real>("K_struct")),
	_tau_0(getParam<Real>("tau_0")),
	
	// Initial values of the state variables
	_init_d_struct(getParam<Real>("init_d_struct")),
    _init_rho_c(getParam<Real>("init_rho_c")),
	_init_rho_w(getParam<Real>("init_rho_w")),
	
	// PSBs variables and parameters
	_f_PSB_0(getParam<Real>("f_PSB_0")),
	_f_PSB_inf(getParam<Real>("f_PSB_inf")),
	_k_PSB(getParam<Real>("k_PSB")),
	_epsilon_p_eff_cum_PSB(getParam<Real>("epsilon_p_eff_cum_PSB")),
	_k_c_PSB(getParam<Real>("k_c_PSB")),
	_y_PSB(getParam<Real>("y_PSB")),
	_f_w_PSB(getParam<Real>("f_w_PSB")),
	_eta_PSB(getParam<Real>("eta_PSB")),
	
	_init_d_struct_PSB(getParam<Real>("init_d_struct_PSB")),
	_init_rho_PSB(getParam<Real>("init_rho_PSB")),
	
	// Interaction matrix coefficients between slip systems
	
	_A_self(getParam<Real>("A_self")),
	_A_copl(getParam<Real>("A_copl")),
	_A_CS(getParam<Real>("A_CS")),
	_A_GJ(getParam<Real>("A_GJ")),
	_A_HL(getParam<Real>("A_HL")),
	_A_LC(getParam<Real>("A_LC")),
	
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
	
    // Wall volume fraction
    _f_w(declareProperty<Real>("f_w")),

	// Channel aspect ratio of dislocation structure
    _eta(declareProperty<Real>("eta")),
    _eta_vector(declareProperty<std::vector<Real>>("eta_vector")),
	
	// Characteristic dislocation structure length
    _d_struct(declareProperty<Real>("d_struct")),
    _d_struct_old(getMaterialPropertyOld<Real>("d_struct")),
	
	// Mean glide distance for dislocations in channel phase
    _l_c(declareProperty<Real>("l_c")),
    _l_c_vector(declareProperty<std::vector<Real>>("l_c_vector")),
    
	// Characteristic PSB length
    _d_struct_PSB(declareProperty<Real>("d_struct_PSB")),
	_d_struct_PSB_old(getMaterialPropertyOld<Real>("d_struct_PSB")),
	
	// Mean glide distance for dislocations in PSB
    _l_PSB(declareProperty<Real>("l_PSB")),
	
	// PSB volume fraction
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

    // Increments of state variables
    _rho_c_increment(_number_slip_systems, 0.0),
	_rho_w_increment(_number_slip_systems, 0.0),
	_rho_PSB_increment(_number_slip_systems, 0.0),
    _backstress_c_increment(_number_slip_systems, 0.0),
	_backstress_w_increment(_number_slip_systems, 0.0),
	_backstress_PSB_increment(_number_slip_systems, 0.0),
	
	// Resize local caching vectors used for substepping
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
  CrystalPlasticityDislocationUpdateBase::initQpStatefulProperties();
  
  Real taylor_hardening;
  
  // Peierls stress for wall and channel
  Real s_0w;
  Real s_0c;

  // Initialize the dislocation density size
  _rho_c[_qp].resize(_number_slip_systems);
  _rho_w[_qp].resize(_number_slip_systems);
  _rho_PSB[_qp].resize(_number_slip_systems);
  
  // Initialize the backstress size
  _backstress_c[_qp].resize(_number_slip_systems);
  _backstress_w[_qp].resize(_number_slip_systems);
  _backstress_PSB[_qp].resize(_number_slip_systems);
  
  // Initialize slip resistance size
  _slip_resistance_c[_qp].resize(_number_slip_systems);
  _slip_resistance_w[_qp].resize(_number_slip_systems);
  _slip_resistance_PSB[_qp].resize(_number_slip_systems);
  
  // Initialize slip increment size
  _slip_increment_c[_qp].resize(_number_slip_systems);
  _slip_increment_w[_qp].resize(_number_slip_systems);
  _slip_increment_PSB[_qp].resize(_number_slip_systems);
  
  // Initialize instantaneous plastic deformation tangent at the slip system level
  _dslip_dtau[_qp].resize(_number_slip_systems);
  
  // Initialize wall volume fraction
  _f_w[_qp] = _f_0;
  
  // Initialize channel aspect ratio of the dislocation structure
  _eta[_qp] = _eta_0;
  _eta_vector[_qp].resize(_number_slip_systems);
  
  if (isParamValid("eta_0_alpha")) {
    for (const auto i : make_range(_number_slip_systems))
    {
      _eta_vector[_qp][i] = _eta_0_alpha[i];
    }
  }
  
  // Initialize characteristic dislocation structure length
  
  if (_read_init_d) { // Read initial characteristic dislocation structure length from file
  
    _d_struct[_qp] = _read_init_d->getData(_current_elem, 0);
    
    } else { // Initialize uniform characteristic dislocation structure length
    
    _d_struct[_qp] = _init_d_struct;
    
  }
  
  // Initialize characteristic PSB length
  _d_struct_PSB[_qp] = _d_struct[_qp];
  
  // Initialize mean glide distance for dislocations in channel phase
  _l_c[_qp] = _eta[_qp] * _d_struct[_qp];
  _l_c_vector[_qp].resize(_number_slip_systems);
  
  if (isParamValid("eta_0_alpha")) {
    for (const auto i : make_range(_number_slip_systems))
    {
      _l_c_vector[_qp][i] = _eta_vector[_qp][i] * _d_struct[_qp];
    }
  }
  
  // Initialize mean glide distance for dislocations in PSB
  _l_PSB[_qp] = _eta_PSB * _d_struct_PSB[_qp];
  
  // Initialize PSB volume fraction
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
  
  if (isParamValid("s_0w") && isParamValid("s_0c")) {
    s_0w = _s_0w;
    s_0c = _s_0c;
  } else {
    s_0w = _s_0;
    s_0c = _s_0;
  }
  
  for (const auto i : make_range(_number_slip_systems))
  {
    // Add Peierls stress
    _slip_resistance_c[_qp][i] = s_0c;
	
	_slip_resistance_c[_qp][i] += (_shear_modulus * _burgers_vector_mag
	                           * std::sqrt(_A_int(i,i) * _rho_c[_qp][i]));
							   
	// Add Peierls stress
	_slip_resistance_w[_qp][i] = s_0w;
	
	taylor_hardening = 0.0;
	
	for (const auto j : make_range(_number_slip_systems))
    {
      taylor_hardening += _A_int(i,j) * _rho_w[_qp][j];
	}
	
    _slip_resistance_w[_qp][i] += _shear_modulus * _burgers_vector_mag
	                           * std::sqrt(taylor_hardening);
	
	_slip_resistance_PSB[_qp][i] = s_0c;
	
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
  Real G0 = _A_self;
  Real G1 = _A_copl;
  Real G2 = _A_CS;
  Real G3 = _A_GJ;
  Real G4 = _A_HL;
  Real G5 = _A_LC;
  
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
        _gamma_o * std::pow(stress_ratio_c, 1.0 / _m_exp);
      
    if (effective_stress_c < 0.0)
      _slip_increment_c[_qp][i] *= -1.0;

  }
  
  // Slip rate inside wall
  
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress_w = _tau[_qp][i] - _backstress_w[_qp][i];
    
    stress_ratio_w = std::abs(effective_stress_w / _slip_resistance_w[_qp][i]);
    
    _slip_increment_w[_qp][i] =
        _gamma_o * std::pow(stress_ratio_w, 1.0 / _m_exp);
      
    if (effective_stress_w < 0.0)
      _slip_increment_w[_qp][i] *= -1.0;

  }
  
  // Slip rate inside PSB
  
  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress_PSB = _tau[_qp][i] - _backstress_PSB[_qp][i];
    
    stress_ratio_PSB = std::abs(effective_stress_PSB / _slip_resistance_PSB[_qp][i]);
    
    _slip_increment_PSB[_qp][i] =
        _gamma_o * std::pow(stress_ratio_PSB, 1.0 / _m_exp);
      
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
  
  // Peierls stress for wall and channel
  Real s_0w;
  Real s_0c;
  
  if (isParamValid("s_0w") && isParamValid("s_0c")) {
    s_0w = _s_0w;
    s_0c = _s_0c;
  } else {
    s_0w = _s_0;
    s_0c = _s_0;
  }
  
  for (const auto i : make_range(_number_slip_systems))
  {
    // Add Peierls stress
    _slip_resistance_c[_qp][i] = s_0c;
	
	_slip_resistance_c[_qp][i] += (_shear_modulus * _burgers_vector_mag
	                           * std::sqrt(_A_int(i,i) * _rho_c[_qp][i]));
							   
	// Add Peierls stress
	_slip_resistance_w[_qp][i] = s_0w;
	
	taylor_hardening = 0.0;
	
	for (const auto j : make_range(_number_slip_systems))
    {
      taylor_hardening += _A_int(i,j) * _rho_w[_qp][j];
	}
	
    _slip_resistance_w[_qp][i] += _shear_modulus * _burgers_vector_mag
	                           * std::sqrt(taylor_hardening);
	
	_slip_resistance_PSB[_qp][i] = s_0c;
	
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

      dslip_dtau_c[i] = _gamma_o / _m_exp *
                      std::pow(stress_ratio_c, 1.0 / _m_exp - 1.0) /
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

      dslip_dtau_w[i] = _gamma_o / _m_exp *
                      std::pow(stress_ratio_w, 1.0 / _m_exp - 1.0) /
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

      dslip_dtau_PSB[i] = _gamma_o / _m_exp *
                      std::pow(stress_ratio_PSB, 1.0 / _m_exp - 1.0) /
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
  
  // Mean glide distance for dislocations in channel phase
  // temporary variable used for each slip system
  Real l_c = _l_c[_qp];

  // channel and wall dislocation density increment
  for (const auto i : make_range(_number_slip_systems))
  {
    if (isParamValid("eta_0_alpha")) {
		
      l_c = _l_c_vector[_qp][i];
	}

    // Multiplication and annihilation
	// note that _slip_increment here is the rate
	// and the rate equation gets multiplied by time step in updateStateVariables
    _rho_c_increment[i] = _k_c / l_c - 2.0 * _y_c * _rho_c[_qp][i];
    _rho_c_increment[i] *= std::abs(_slip_increment_c[_qp][i]) / _burgers_vector_mag;
	
	slip_increment_matrix = (1.0 - _f_w[_qp]) * _slip_increment_c[_qp][i] + _f_w[_qp] * _slip_increment_w[_qp][i];
	
	_rho_w_increment[i] = ( _k_w * std::abs(_slip_increment_c[_qp][i]) ) / ( _burgers_vector_mag * l_c );
    _rho_w_increment[i] -= ( 2.0 * _y_w * _rho_w[_qp][i] * std::abs(slip_increment_matrix) ) / _burgers_vector_mag;
	
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
  
  // secant elastoplastic Poisson's ratio
  Real nu_p;

  // S1212 is a component of the Eshelby tensor
  Real S_1212;
  Real S_1212_PSB;
  
  // secant elastoplastic shear modulus
  Real mu_m;
	
  // Function that depends on the channel aspect ratio
  Real f_eta;
  Real f_eta_PSB;
  
  // channel aspect ratio 
  // temporary variable for each slip system
  Real eta = _eta[_qp];
  
  K_shape = eta * std::sqrt( std::pow(eta, 2.0) - 1.0) - std::acosh(eta);
  K_shape *= 2.0 * libMesh::pi * eta;
  K_shape /= std::sqrt(std::pow(std::pow(eta, 2.0) - 1.0, 3.0));
  
  K_shape_PSB = _eta_PSB * std::sqrt( std::pow(_eta_PSB, 2.0) - 1.0) - std::acosh(_eta_PSB);
  K_shape_PSB *= 2.0 * libMesh::pi * _eta_PSB;
  K_shape_PSB /= std::sqrt(std::pow(std::pow(_eta_PSB, 2.0) - 1.0, 3.0));
  
  for (const auto i : make_range(_number_slip_systems)) 
  {
    f_accom = _dslip_dtau[_qp][i] * _substep_dt / 2.0;
    
    nu_p = _nu + (2.0/3.0) * (1.0 + _nu) * _shear_modulus * f_accom;
	nu_p /= 1.0 + (4.0/3.0) * (1.0 + _nu) * _shear_modulus * f_accom;
	
	if (isParamValid("eta_0_alpha")) {
		
	  eta = _eta_vector[_qp][i];

      K_shape = eta * std::sqrt( std::pow(eta, 2.0) - 1.0) - std::acosh(eta);
      K_shape *= 2.0 * libMesh::pi * eta;
      K_shape /= std::sqrt(std::pow(std::pow(eta, 2.0) - 1.0, 3.0));  
	}
	
	S_1212 = std::pow(eta, 2.0) - 1.75 - 2.0 * nu_p * std::pow(eta, 2.0) + nu_p;
	S_1212 *= K_shape;
	S_1212 += libMesh::pi * std::pow(eta, 2.0);
	S_1212 /= 8.0 * libMesh::pi * (1.0 - nu_p) * (std::pow(eta, 2.0) - 1.0);
	
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

// Update dislocation densities and backstress based on increments
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

// Calculate channel aspect ratio based on (41), (42), (43)
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

// Calculate characteristic dislocation structure length based on (40)
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

// Calculate PSB volume fraction
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
