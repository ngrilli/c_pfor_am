// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// University of Bristol
// 04 November 2025

#include "CrystalPlasticityDislocationUpdateNi.h"
#include "libmesh/int_range.h"
#include <cmath>
#include "Function.h"

registerMooseObject("c_pfor_amApp", CrystalPlasticityDislocationUpdateNi);

InputParameters
CrystalPlasticityDislocationUpdateNi::validParams()
{
  InputParameters params = CrystalPlasticityDislocationUpdateBase::validParams();
  params.addClassDescription("Dislocation based model for crystal plasticity "
                             "using the stress update code for Ni-based alloy (IN718). "
                             "Includes slip, creep and backstress, solid solution strenthening and precipitation hardening. ");
  params.addParam<Real>("ao", 0.001, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");  
  params.addParam<MaterialPropertyName>("xm_matprop",
    "Optional xm material property for exponent for slip rate. ");
  params.addParam<FunctionName>("ao_function",
    "Optional function for slip prefactor. If provided, the slip prefactor can be set as a function of time. "
    "This is useful for an initial plastic deformation followed by creep load. ");
  params.addParam<bool>("use_kocks_T_dependence_for_xm", false, "Use Kocks 1976 temperature dependence for xm. ");
  params.addParam<bool>("creep_activated", false, "Activate creep strain rate.");
  params.addParam<Real>("creep_ao", 0.0, "creep rate coefficient");
  params.addParam<Real>("creep_xm", 0.1, "exponent for creep rate");
  params.addParam<FunctionName>("creep_ao_function",
    "Optional function for creep prefactor. If provided, the creep prefactor can be set as a function of time. "
    "This is useful for an initial plastic deformation followed by creep load. ");
  params.addParam<FunctionName>("creep_resistance_function",
    "Optional function for creep resistance. If provided, the creep resistance can be set as a function of time. "
    "This is useful for differentiating resistance for slip and creep. ");
  params.addParam<Real>("m_exponent", 0.0, "Exponent on time in power-law equation");
  params.addParam<Real>("creep_t0", 0.0, "Initial time for tertiary creep");
  params.addParam<Real>("creep_t_denominator", 1.0, "Denominator for the tertiary creep law");
  params.addParam<bool>("cap_slip_increment", false, "Cap the absolute value of the slip increment "
                                                     "in one time step to _slip_incr_tol. ");
  params.addParam<Real>("burgers_vector_mag",0.000256,"Magnitude of the Burgers vector");
  params.addParam<Real>("shear_modulus",86000.0,"Shear modulus in Taylor hardening law G");
  params.addParam<Real>("dshear_dT", 0.0,"Temperature dependence of shear modulus (C44)");  // NEW
  params.addParam<Real>("alpha_0",0.3,"Prefactor of Taylor hardening law, alpha");
  //params.addParam<Real>("r", 1.4, "Latent hardening coefficient");
  params.addParam<Real>("r_self", 1.0, "Self-hardening coefficient (i=j)"); // NEW
  params.addParam<Real>("r_oct_oct", 1.4, "Latent hardening {111} <-> {111} (e.g., 1.4)"); // NEW
  params.addParam<Real>("r_cub_cub", 1.0, "Latent hardening {100} <-> {100} (e.g., 1.0)"); // NEW
  params.addParam<Real>("r_oct_cub", 1.2, "Latent hardening {111} <-> {100} (e.g., 1.2)"); // NEW  
  params.addParam<Real>("tau_c_0_oct", 0.112, "Peierls stress for {111} octahedral slip");
  params.addParam<Real>("tau_c_0_cub", 1700.0, "Peierls stress for {100} cubic slip (set high to disable at RT)");
  params.addParam<FunctionName>("tau_c_0_function",
    "Optional function for Peierls stress. If provided, the Peierls stress can be set as a function of time. "
    "This is useful for time dependent solid solution strenthening and precipitation hardening, e.g. during thermal ageing. ");
  params.addParam<Real>("tau_ss", 0.0, "Solid solution strengthening contribution"); // <--NEW
  params.addParam<Real>("k_hp", 0.0, "Hall-Petch coefficient (e.g., in MPa*sqrt(um))"); // <--NEW
  // --- Parameters for Precipitate Strengthening (PDF 3.2) --- 
  params.addParam<Real>("C_g_prime", 0.0, "C_gamma' constant in strengthening eqn"); // NEW
  params.addParam<Real>("Gamma_APB_g_prime", 0.0, "APB energy for g' (L12) [J/m^2]"); // NEW
  params.addParam<Real>("f_vol_g_prime", 0.0, "Volume fraction of g' (L12)"); // NEW
  params.addParam<Real>("C_g_pp", 0.0, "C_gamma'' constant in strengthening eqn"); // NEW
  params.addParam<Real>("Gamma_APB_g_pp", 0.0, "APB energy for g'' (DO22) [J/m^2]"); // NEW
  params.addParam<Real>("f_vol_g_pp", 0.0, "Volume fraction of g'' (DO22)"); // NEW
  params.addParam<Real>("init_r_eff_g_prime", 0.0, "Initial effective radius of g' (L12) [m]"); // NEW
  params.addParam<Real>("init_r_eff_g_pp", 0.0, "Initial effective radius of g'' (DO22) [m]"); // NEW
  params.addParam<Real>("C_shear_g_prime", 0.0, "Shearing efficiency coefficient for g' softening"); // NEW
  params.addParam<Real>("C_shear_g_pp", 0.0, "Shearing efficiency coefficient for g'' softening");  // NEW
  //Parameters for DRV at elevated temperatures  
  params.addParam<Real>("k_0",100.0,"Coefficient K in SSD evolution, representing accumulation rate");
  params.addParam<Real>("y_c",0.0026,"Critical annihilation diameter");
  params.addParam<Real>("Q_drv", 150.0e3, "Activation energy for dynamic recovery (e.g., J/mol)"); // NEW
  params.addParam<Real>("R_gas_constant", 8.314, "Gas constant (e.g., J/mol/K)");  // NEW
  params.addParam<Real>("h",0.0,"Direct hardening coefficient for backstress");
  params.addParam<Real>("h_D",0.0,"Dynamic recovery coefficient for backstress");
 // --- NEW: Parameters for Intragranular Backstress (Agaram et al.) ---
  params.addParam<Real>("k_52", 0.0, "Slope parameter for intragranular backstress"); // NEW
  params.addParam<Real>("k_32", 0.0, "Saturation parameter for intragranular backstress"); // NEW
  params.addParam<Real>("k_D", 0.0, "Precipitate dislocation generation parameter (for Z_o)"); // NEW
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
  params.addParam<UserObjectName>("read_grain_size",
                                  "The ElementPropertyReadFile "
                                  "GeneralUserObject to read element value "
                                  "of the grain size (scalar)");
  params.addCoupledVar("dslip_increment_dedge",0.0,"Directional derivative of the slip rate along the edge motion direction.");
  params.addCoupledVar("dslip_increment_dscrew",0.0,"Directional derivative of the slip rate along the screw motion direction.");
  params.addCoupledVar("temperature", 303.0,"Temperature, initialize at room temperature");
  params.addParam<Real>("reference_temperature",303.0,"reference temperature for thermal expansion");
  params.addParam<Real>("dCRSS_dT_A",1.0,"A coefficient for the exponential decrease of the CRSS of octahedral slip"
                        "with temperature: A + B exp(- C * (T - 303.0))");
  params.addParam<Real>("dCRSS_dT_B",0.0,"B coefficient for the exponential decrease of the CRSS of octahedral slip "
                        "with temperature: A + B exp(- C * (T - 303.0))");
  params.addParam<Real>("dCRSS_dT_C",0.0,"C coefficient for the exponential decrease of the CRSS of octahedral slip"
                        "with temperature: A + B exp(- C * (T - 303.0))");
  params.addParam<Real>("dCRSS_dT_A_cub", 1.0, "A coeff for Cubic {100} slip T-dependence"); // NEW
  params.addParam<Real>("dCRSS_dT_B_cub", 0.0, "B coeff for Cubic {100} slip T-dependence"); // NEW
  params.addParam<Real>("dCRSS_dT_C_cub", 0.0, "C coeff for Cubic {100} slip T-dependence"); // NEW   
  params.addParam<Real>("r_prep_tol", 1.0e-6, "Tolerance of precipitation radii evolution"); // NEW                          
  return params;
}

CrystalPlasticityDislocationUpdateNi::CrystalPlasticityDislocationUpdateNi(
    const InputParameters & parameters)
  : CrystalPlasticityDislocationUpdateBase(parameters),
  
    // Constitutive model parameters
    _ao(getParam<Real>("ao")),
    _xm(getParam<Real>("xm")),
    _include_xm_matprop(parameters.isParamValid("xm_matprop")),
    _xm_matprop(_include_xm_matprop
                ? &getMaterialProperty<Real>("xm_matprop")
                : nullptr),
    _ao_function(this->isParamValid("ao_function")
                ? &this->getFunction("ao_function")
                : NULL),
    _use_kocks_T_dependence_for_xm(getParam<bool>("use_kocks_T_dependence_for_xm")),
    _creep_activated(getParam<bool>("creep_activated")),
    _creep_ao(getParam<Real>("creep_ao")),
    _creep_xm(getParam<Real>("creep_xm")),
    _creep_ao_function(this->isParamValid("creep_ao_function")
                       ? &this->getFunction("creep_ao_function")
                       : NULL),
    _creep_resistance_function(this->isParamValid("creep_resistance_function")
                       ? &this->getFunction("creep_resistance_function")
                       : NULL),
    _m_exponent(getParam<Real>("m_exponent")),
    _creep_t0(getParam<Real>("creep_t0")),
    _creep_t_denominator(getParam<Real>("creep_t_denominator")),
    _cap_slip_increment(getParam<bool>("cap_slip_increment")),
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")),
	_shear_modulus(getParam<Real>("shear_modulus")),
  _dshear_dT(getParam<Real>("dshear_dT")),         // NEW 
	_alpha_0(getParam<Real>("alpha_0")),
   // _r(getParam<Real>("r")),
    _r_self(getParam<Real>("r_self")),       // NEW
    _r_oct_oct(getParam<Real>("r_oct_oct")), // NEW
    _r_cub_cub(getParam<Real>("r_cub_cub")), // NEW
    _r_oct_cub(getParam<Real>("r_oct_cub")), // NEW 
  _tau_c_0_oct(getParam<Real>("tau_c_0_oct")), // NEW
	_tau_c_0_cub(getParam<Real>("tau_c_0_cub")), // NEW
	_tau_c_0_function(this->isParamValid("tau_c_0_function")
                       ? &this->getFunction("tau_c_0_function")
                       : NULL),
  _tau_ss(getParam<Real>("tau_ss")), // new 
  _k_hp(getParam<Real>("k_hp")),       // new 
    // --- Init Parameters for Precipitate Strengthening --- 
    _C_g_prime(getParam<Real>("C_g_prime")),               // NEW
    _Gamma_APB_g_prime(getParam<Real>("Gamma_APB_g_prime")), // NEW
    _f_vol_g_prime(getParam<Real>("f_vol_g_prime")),       // NEW
    _C_g_pp(getParam<Real>("C_g_pp")),                      // NEW
    _Gamma_APB_g_pp(getParam<Real>("Gamma_APB_g_pp")),      // NEW
    _f_vol_g_pp(getParam<Real>("f_vol_g_pp")),              // NEW

	_k_0(getParam<Real>("k_0")),
	_y_c(getParam<Real>("y_c")),
  _Q_drv(getParam<Real>("Q_drv")), // new
  _R_gas_constant(getParam<Real>("R_gas_constant")),	// new
  
    _init_r_eff_g_prime(getParam<Real>("init_r_eff_g_prime")), // NEW
    _init_r_eff_g_pp(getParam<Real>("init_r_eff_g_pp")),      // NEW   
    _C_shear_g_prime(getParam<Real>("C_shear_g_prime")),      // NEW 
    _C_shear_g_pp(getParam<Real>("C_shear_g_pp")),         // NEW   
  
	// Backstress parameters
	_h(getParam<Real>("h")),
	_h_D(getParam<Real>("h_D")),
  // --- NEW: Initialize Intragranular Backstress Parameters 
  _k_52(getParam<Real>("k_52")),
  _k_32(getParam<Real>("k_32")),
  _k_D(getParam<Real>("k_D")),	
	// Initial values of the state variables
    _init_rho_ssd(getParam<Real>("init_rho_ssd")),
    _init_rho_gnd_edge(getParam<Real>("init_rho_gnd_edge")),
    _init_rho_gnd_screw(getParam<Real>("init_rho_gnd_screw")),
	
	// Tolerance on dislocation density update
	_rho_tol(getParam<Real>("rho_tol")),
	
	// State variables of the dislocation model
	
    _rho_ssd(declareProperty<std::vector<Real>>("rho_ssd")),
    _rho_ssd_old(getMaterialPropertyOld<std::vector<Real>>("rho_ssd")),
    _rho_gnd_edge(declareProperty<std::vector<Real>>("rho_gnd_edge")),
   	_rho_gnd_edge_old(getMaterialPropertyOld<std::vector<Real>>("rho_gnd_edge")),
  	_rho_gnd_screw(declareProperty<std::vector<Real>>("rho_gnd_screw")),
    _rho_gnd_screw_old(getMaterialPropertyOld<std::vector<Real>>("rho_gnd_screw")),
  
    // --- Declare New ISVs for Precipitate Radius --- 
    _r_eff_g_prime(declareProperty<std::vector<Real>>("r_eff_g_prime")), // NEW
    _r_eff_g_prime_old(getMaterialPropertyOld<std::vector<Real>>("r_eff_g_prime")), // NEW
    _r_eff_g_pp(declareProperty<std::vector<Real>>("r_eff_g_pp")),    // NEW
    _r_eff_g_pp_old(getMaterialPropertyOld<std::vector<Real>>("r_eff_g_pp")),    // NEW  
    // Backstress variable
    _backstress_inter(declareProperty<std::vector<Real>>("backstress_inter")), // RENAMED
    _backstress_inter_old(getMaterialPropertyOld<std::vector<Real>>("backstress_inter")), // RENAMED
    _backstress_intra(declareProperty<std::vector<Real>>("backstress_intra")), // NEW
    _backstress_intra_old(getMaterialPropertyOld<std::vector<Real>>("backstress_intra")), // NEW

    // increments of state variables
	
    _rho_ssd_increment(_number_slip_systems, 0.0),
    _rho_gnd_edge_increment(_number_slip_systems, 0.0),
    _rho_gnd_screw_increment(_number_slip_systems, 0.0),
    _backstress_inter_increment(_number_slip_systems, 0.0), // RENAMED
    _backstress_intra_increment(_number_slip_systems, 0.0), // NEW
	  
    _r_eff_g_prime_increment(_number_slip_systems, 0.0), // NEW
    _r_eff_g_pp_increment(_number_slip_systems, 0.0),    // NEW

	// resize local caching vectors used for substepping
    _previous_substep_rho_ssd(_number_slip_systems, 0.0),
	 _previous_substep_rho_gnd_edge(_number_slip_systems, 0.0),
	 _previous_substep_rho_gnd_screw(_number_slip_systems, 0.0),
   _previous_substep_backstress_inter(_number_slip_systems, 0.0), // RENAMED
	 _previous_substep_backstress_intra(_number_slip_systems, 0.0), // NEW
    _previous_substep_r_eff_g_prime(_number_slip_systems, 0.0), // NEW
    _previous_substep_r_eff_g_pp(_number_slip_systems, 0.0),    // NEW

    _rho_ssd_before_update(_number_slip_systems, 0.0),
    _rho_gnd_edge_before_update(_number_slip_systems, 0.0),
    _rho_gnd_screw_before_update(_number_slip_systems, 0.0),  	
    _backstress_inter_before_update(_number_slip_systems, 0.0), // RENAMED
    _backstress_intra_before_update(_number_slip_systems, 0.0), // NEW

    _r_eff_g_prime_before_update(_number_slip_systems, 0.0), // NEW
    _r_eff_g_pp_before_update(_number_slip_systems, 0.0),    // NEW

    // Twinning contributions, if used
    _include_twinning_in_Lp(parameters.isParamValid("total_twin_volume_fraction")),
    _twin_volume_fraction_total(_include_twinning_in_Lp
                                    ? &getMaterialPropertyOld<Real>("total_twin_volume_fraction")
                                    : nullptr),

    // UserObject to read the initial GND density from file						
    _read_initial_gnd_density(isParamValid("read_initial_gnd_density")
                               ? &getUserObject<ElementPropertyReadFile>("read_initial_gnd_density")
                               : nullptr),

    // UserObject to read the grain size from EBSD file	                          
		_read_grain_size(isParamValid("read_grain_size")
                               ? &getUserObject<ElementPropertyReadFile>("read_grain_size")  // NEW
                               : nullptr),				

    // Directional derivatives of the slip rate
    _include_slip_gradients(isParamValid("dslip_increment_dedge") && isParamValid("dslip_increment_dscrew")),
    _dslip_increment_dedge(_include_slip_gradients
                            ? coupledArrayValue("dslip_increment_dedge")
                            : _default_array_value_zero),
    _dslip_increment_dscrew(_include_slip_gradients
                            ? coupledArrayValue("dslip_increment_dscrew")
                            : _default_array_value_zero),
	
	// Temperature dependent properties
	_temperature(coupledValue("temperature")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _dCRSS_dT_A(getParam<Real>("dCRSS_dT_A")),
	_dCRSS_dT_B(getParam<Real>("dCRSS_dT_B")),
	_dCRSS_dT_C(getParam<Real>("dCRSS_dT_C")),
    _dCRSS_dT_A_cub(getParam<Real>("dCRSS_dT_A_cub")), // NEW
    _dCRSS_dT_B_cub(getParam<Real>("dCRSS_dT_B_cub")), // NEW
    _dCRSS_dT_C_cub(getParam<Real>("dCRSS_dT_C_cub")),  // NEW
    _r_prep_tol(getParam<Real>("r_prep_tol")),  // NEW
    // store edge and screw slip directions to calculate directional derivatives
    // of the plastic slip rate	
    _edge_slip_direction(declareProperty<std::vector<Real>>("edge_slip_direction")),
	 _screw_slip_direction(declareProperty<std::vector<Real>>("screw_slip_direction"))
{
}

void
CrystalPlasticityDislocationUpdateNi::initQpStatefulProperties()
{
  // Slip resistance is resized here
  CrystalPlasticityDislocationUpdateBase::initQpStatefulProperties();
  
  Real taylor_hardening;
  
  // Temperature dependence of the CRSS
  Real temperature_dependence;

  // Initialize the dislocation density size
  _rho_ssd[_qp].resize(_number_slip_systems);
  _rho_gnd_edge[_qp].resize(_number_slip_systems);
  _rho_gnd_screw[_qp].resize(_number_slip_systems);
  
  // Initialize the backstress size
  _backstress_inter[_qp].resize(_number_slip_systems); // RENAMED
  _backstress_intra[_qp].resize(_number_slip_systems); // NEW
  _r_eff_g_prime[_qp].resize(_number_slip_systems); // NEW
  _r_eff_g_pp[_qp].resize(_number_slip_systems);    // NEW

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
	
	_backstress_inter[_qp][i] = 0.0; // RENAMED
	_backstress_intra[_qp][i] = 0.0; // NEW

     // Init precipitate radii ISVs 
    _r_eff_g_prime[_qp][i] = _init_r_eff_g_prime; // NEW
    _r_eff_g_pp[_qp][i] = _init_r_eff_g_pp;       // NEW 
  }

  Real tau_hp = 0.0;
  if (_read_grain_size)
  {
      Real d_elem = _read_grain_size->getData(_current_elem, 0); // NEW
      if (d_elem > 1.0e-9) // avoid zero
         tau_hp = _k_hp / std::sqrt(d_elem);
  }

  // Calculate temperature dependence for all systems
  for (const auto i : make_range(_number_slip_systems))
  {
    Real tau_c_0_i; // NEW: System-dependent Peierls stress
    // We assume first 12 systems are Octahedral {111}
    if (i < 12) // NEW
    { tau_c_0_i = _tau_c_0_oct; // NEW
  // Critical resolved shear stress decreases exponentially with temperature
  // A + B exp(- C * (T - T0)) 
       temperature_dependence = ( _dCRSS_dT_A + _dCRSS_dT_B 
                                             * std::exp(- _dCRSS_dT_C * (_temperature[_qp] - _reference_temperature)));
    } // NEW
    // We assume systems 12-17 are Cubic {100}
    else // NEW
    { tau_c_0_i = _tau_c_0_cub; // NEW
        temperature_dependence = ( _dCRSS_dT_A_cub + _dCRSS_dT_B_cub // NEW
                                 * std::exp(- _dCRSS_dT_C_cub * (_temperature[_qp] - _reference_temperature))); // NEW
    } // NEW

    // Ensure non-zero/negative values to prevent numerical failure
    if (temperature_dependence <= 0) // NEW
        temperature_dependence = 1.0e-12; // NEW 
 

  // Initialize value of the slip resistance
  // as a function of the dislocation density

    // Add Peierls stress AND solid solution strengthening
    _slip_resistance[_qp][i] = tau_c_0_i + _tau_ss;

    taylor_hardening = 0.0;
	  
    for (const auto j : make_range(_number_slip_systems))
    {
      // Use new anisotropic latent hardening logic
      // This replaces the simple `iplane == jplane` logic

      Real r_ij = 0.0; // NEW: Interaction coefficient
      
      if (i == j) // 1. Self-hardening
      { // NEW
          r_ij = _r_self; // NEW
      } // NEW
      else if (i < 12 && j < 12) // 2. Octahedral <-> Octahedral
      { // NEW
          r_ij = _r_oct_oct; // NEW
      } // NEW
      else if (i >= 12 && j >= 12) // 3. Cubic <-> Cubic
      { // NEW
          r_ij = _r_cub_cub; // NEW
      } // NEW
      else // 4. Octahedral <-> Cubic
      { // NEW
          r_ij = _r_oct_cub; // NEW
      } // NEW

      taylor_hardening += r_ij * (_rho_ssd[_qp][j]  // NEW
                          + std::abs(_rho_gnd_edge[_qp][j]) // NEW
                          + std::abs(_rho_gnd_screw[_qp][j])); // NEW
	  }
  
	  Real G_T = _shear_modulus * (1.0 - _dshear_dT * (_temperature[_qp] - _reference_temperature)); // NEW  
    Real tau_dislo = _alpha_0 * G_T * _burgers_vector_mag * std::sqrt(taylor_hardening); // NEW: dislocation part

    // Calculate Precipitate Strengthening
    // Saeede Ghorbanpour, M Irene J. Beyerlein, Marko Knezevic, et al.
    //A crystal plasticity model incorporating the effects of precipitates in superalloys: Application to tensile, compressive, and cyclic deformation of Inconel 718,
    // International Journal of Plasticity, 2017

    Real r_g_prime_i = _r_eff_g_prime[_qp][i]; // NEW
    Real r_g_pp_i = _r_eff_g_pp[_qp][i];       // NEW
    Real tau_g_prime = 0.0; // NEW
    if (_f_vol_g_prime > 0) // NEW
      tau_g_prime = _C_g_prime * std::pow(_Gamma_APB_g_prime / _burgers_vector_mag, 1.5) * std::pow(_f_vol_g_prime * r_g_prime_i /_burgers_vector_mag, 0.5); // NEW
    
    Real tau_g_pp = 0.0; // NEW
    if (_f_vol_g_pp > 0) // NEW
      tau_g_pp = _C_g_pp * std::pow(_Gamma_APB_g_pp / _burgers_vector_mag, 1.5) * std::pow(_f_vol_g_pp * r_g_pp_i / _burgers_vector_mag, 0.5); // NEW

    Real tau_precip_i = tau_g_prime + tau_g_pp; // NEW
    // Add thermally-dependent strengthening components
	//_slip_resistance[_qp][i] += (_alpha_0 * _shear_modulus * _burgers_vector_mag * std::sqrt(taylor_hardening) * temperature_dependence);
	  _slip_resistance[_qp][i] += tau_precip_i + (tau_hp + tau_dislo) * temperature_dependence ;	 // NEW
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
CrystalPlasticityDislocationUpdateNi::calculateSchmidTensor(
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
CrystalPlasticityDislocationUpdateNi::setInitialConstitutiveVariableValues()
{
  // Initialize state variables with the value at the previous time step
  _rho_ssd[_qp] = _rho_ssd_old[_qp];
  _previous_substep_rho_ssd = _rho_ssd_old[_qp];
  _rho_gnd_edge[_qp] = _rho_gnd_edge_old[_qp];
  _previous_substep_rho_gnd_edge = _rho_gnd_edge_old[_qp];
  _rho_gnd_screw[_qp] = _rho_gnd_screw_old[_qp];
  _previous_substep_rho_gnd_screw = _rho_gnd_screw_old[_qp];
  _backstress_inter[_qp] = _backstress_inter_old[_qp]; // RENAMED
  _previous_substep_backstress_inter = _backstress_inter_old[_qp]; // RENAMED
  _backstress_intra[_qp] = _backstress_intra_old[_qp]; // NEW
  _previous_substep_backstress_intra = _backstress_intra_old[_qp]; // NEW
  _r_eff_g_prime[_qp] = _r_eff_g_prime_old[_qp]; // NEW
  _previous_substep_r_eff_g_prime = _r_eff_g_prime_old[_qp]; // NEW
  _r_eff_g_pp[_qp] = _r_eff_g_pp_old[_qp]; // NEW
  _previous_substep_r_eff_g_pp = _r_eff_g_pp_old[_qp]; // NEW  
}

void
CrystalPlasticityDislocationUpdateNi::setSubstepConstitutiveVariableValues()
{
  // Inialize state variable of the next substep
  // with the value at the previous substep
  _rho_ssd[_qp] = _previous_substep_rho_ssd;
  _rho_gnd_edge[_qp] = _previous_substep_rho_gnd_edge;
  _rho_gnd_screw[_qp] = _previous_substep_rho_gnd_screw;
  _backstress_inter[_qp] = _previous_substep_backstress_inter; // RENAMED
  _backstress_intra[_qp] = _previous_substep_backstress_intra; // NEW
  _r_eff_g_prime[_qp] = _previous_substep_r_eff_g_prime; // NEW
  _r_eff_g_pp[_qp] = _previous_substep_r_eff_g_pp;       // NEW  
}

// Slip resistance can be calculated from dislocation density here only
// because it is the first method in which it is used,
// while calculateConstitutiveSlipDerivative is called afterwards
bool
CrystalPlasticityDislocationUpdateNi::calculateSlipRate()
{
  calculateSlipResistance();
  
  // Ratio between effective stress and CRSS
  // temporary variable for each slip system
  Real stress_ratio;
  
  // Difference between RSS and backstress
  // temporary variable for each slip system
  Real effective_stress;
  
  // Slip prefactor: if function is not given
  // the constant value is used
  Real ao;
  
  // Creep prefactor: if function is not given
  // the constant value is used
  Real creep_ao;
  
  // Tertiary creep contribution
  Real tertiary_creep = 1.0;
  
  if (_t > _creep_t0 && _creep_activated) {
    tertiary_creep += std::pow((_t - _creep_t0)/_creep_t_denominator, _m_exponent);
  }
  
  if (_ao_function) {
	  
    ao = _ao_function->value(_t, _q_point[_qp]);
	  
  } else {
	  
    ao = _ao;
	  
  }  
  
  if (_creep_ao_function) {
	  
    creep_ao = _creep_ao_function->value(_t, _q_point[_qp]);
	  
  } else {
	  
    creep_ao = _creep_ao;
	  
  }
  
  // Strain rate sensitivity: if material property is not given
  // the constant value is used
  Real xm;
  
  if (_include_xm_matprop) {
	  
    xm = (*_xm_matprop)[_qp];	  
	  
  } else {
	  
    xm = _xm;
	  
  }
  
  for (const auto i : make_range(_number_slip_systems))
  {
     effective_stress = _tau[_qp][i] - (_backstress_inter[_qp][i] + _backstress_intra[_qp][i]);
    
    stress_ratio = std::abs(effective_stress / _slip_resistance[_qp][i]);
    
    _slip_increment[_qp][i] = ao * std::pow(stress_ratio, 1.0 / xm);
        
    if (_creep_activated) { // add creep rate
      
      if (_creep_resistance_function) {
        stress_ratio = std::abs(effective_stress / _creep_resistance_function->value(_t, _q_point[_qp]));
      }
      
      _slip_increment[_qp][i] += creep_ao * std::pow(stress_ratio, 1.0 / _creep_xm) * tertiary_creep;
    }
      
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
CrystalPlasticityDislocationUpdateNi::calculateSlipResistance()
{
  Real taylor_hardening;
  
  // Temperature dependence of the CRSS
  Real temperature_dependence;
  
  Real tau_hp = 0.0;
  if (_read_grain_size)
  {
      Real d_elem = _read_grain_size->getData(_current_elem, 0); // NEW
      if (d_elem > 1.0e-9) // avoid zero
         tau_hp = _k_hp / std::sqrt(d_elem);
  }  
  // Calculate temperature dependence for all systems
  for (const auto i : make_range(_number_slip_systems))
  {
    Real tau_c_0_i; // NEW: System-dependent Peierls stress
    // We assume first 12 systems are Octahedral {111}
    if (i < 12) // NEW
    { tau_c_0_i = _tau_c_0_oct; // NEW
  // Critical resolved shear stress decreases exponentially with temperature
  // A + B exp(- C * (T - T0)) 
       temperature_dependence = ( _dCRSS_dT_A + _dCRSS_dT_B 
                                             * std::exp(- _dCRSS_dT_C * (_temperature[_qp] - _reference_temperature)));
    } // NEW
    // We assume systems 12-17 are Cubic {100}
    else // NEW
    { tau_c_0_i = _tau_c_0_cub; // NEW
        temperature_dependence = ( _dCRSS_dT_A_cub + _dCRSS_dT_B_cub // NEW
                                 * std::exp(- _dCRSS_dT_C_cub * (_temperature[_qp] - _reference_temperature))); // NEW
    } // NEW

    // Ensure non-zero/negative values to prevent numerical failure
    if (temperature_dependence <= 0) // NEW
        temperature_dependence = 1.0e-12; // NEW 

    // Add Peierls stress
    _slip_resistance[_qp][i] = tau_c_0_i + _tau_ss;

    taylor_hardening = 0.0;
	  
    for (const auto j : make_range(_number_slip_systems))
    {
      // Use new anisotropic latent hardening logic
      // This replaces the simple `iplane == jplane` logic

      Real r_ij = 0.0; // NEW: Interaction coefficient
      
      if (i == j) // 1. Self-hardening
      { // NEW
          r_ij = _r_self; // NEW
      } // NEW
      else if (i < 12 && j < 12) // 2. Octahedral <-> Octahedral
      { // NEW
          r_ij = _r_oct_oct; // NEW
      } // NEW
      else if (i >= 12 && j >= 12) // 3. Cubic <-> Cubic
      { // NEW
          r_ij = _r_cub_cub; // NEW
      } // NEW
      else // 4. Octahedral <-> Cubic
      { // NEW
          r_ij = _r_oct_cub; // NEW
      } // NEW

      taylor_hardening += r_ij * (_rho_ssd[_qp][j]  // NEW
                          + std::abs(_rho_gnd_edge[_qp][j]) // NEW
                          + std::abs(_rho_gnd_screw[_qp][j])); // NEW
	  }
    
		Real G_T = _shear_modulus * (1.0 - _dshear_dT * (_temperature[_qp] - _reference_temperature)); // NEW  
    Real tau_dislo = _alpha_0 * G_T * _burgers_vector_mag * std::sqrt(taylor_hardening); // NEW: dislocation part

    Real r_g_prime_i = _r_eff_g_prime[_qp][i]; // NEW
    Real r_g_pp_i = _r_eff_g_pp[_qp][i];       // NEW
    Real tau_g_prime = 0.0; // NEW
    if (_f_vol_g_prime > 0) // NEW
      tau_g_prime = _C_g_prime * std::pow(_Gamma_APB_g_prime / _burgers_vector_mag, 1.5) * std::pow(_f_vol_g_prime * r_g_prime_i / _burgers_vector_mag, 0.5); // NEW
    
    Real tau_g_pp = 0.0; // NEW
    if (_f_vol_g_pp > 0) // NEW
      tau_g_pp = _C_g_pp * std::pow(_Gamma_APB_g_pp / _burgers_vector_mag, 1.5) * std::pow(_f_vol_g_pp * r_g_pp_i / _burgers_vector_mag, 0.5); // NEW

    Real tau_precip_i = tau_g_prime + tau_g_pp;
    // Add thermally-dependent strengthening components
	//_slip_resistance[_qp][i] += (_alpha_0 * _shear_modulus * _burgers_vector_mag * std::sqrt(taylor_hardening) * temperature_dependence);
	  _slip_resistance[_qp][i] += tau_precip_i + (tau_hp + tau_dislo) * temperature_dependence;	 // NEW
	
  }

}

void
CrystalPlasticityDislocationUpdateNi::calculateEquivalentSlipIncrement(
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
CrystalPlasticityDislocationUpdateNi::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  // Ratio between effective stress and CRSS
  // temporary variable for each slip system
  Real stress_ratio;
  
  // Difference between RSS and backstress
  // temporary variable for each slip system
  Real effective_stress;
  
  // Slip prefactor: if function is not given
  // the constant value is used
  Real ao;
  
  // Creep prefactor: if function is not given
  // the constant value is used
  Real creep_ao;
  
  // Tertiary creep contribution
  Real tertiary_creep = 1.0;
  
  // Creep resistance
  Real creep_resistance;
  
  if (_t > _creep_t0 && _creep_activated) {
    tertiary_creep += std::pow((_t - _creep_t0)/_creep_t_denominator, _m_exponent);
  }
  
  if (_ao_function) {
	  
    ao = _ao_function->value(_t, _q_point[_qp]);
	  
  } else {
	  
    ao = _ao;
	  
  }  
  
  if (_creep_ao_function) {
	  
    creep_ao = _creep_ao_function->value(_t, _q_point[_qp]);
	  
  } else {
	  
    creep_ao = _creep_ao;
	  
  }	
  
  // Strain rate sensitivity: if material property is not given
  // the constant value is used
  Real xm;
  
  if (_include_xm_matprop) {
	  
    xm = (*_xm_matprop)[_qp];	  
	  
  } else {
	  
    xm = _xm;
	  
  }
	
  for (const auto i : make_range(_number_slip_systems))
  {
     effective_stress = _tau[_qp][i] - (_backstress_inter[_qp][i] + _backstress_intra[_qp][i]);
	  
    if (MooseUtils::absoluteFuzzyEqual(effective_stress, 0.0)) {
		
      dslip_dtau[i] = 0.0;
      		
	} else {
		
	  stress_ratio = std::abs(effective_stress / _slip_resistance[_qp][i]);

      dslip_dtau[i] = ao / xm *
                      std::pow(stress_ratio, 1.0 / xm - 1.0) /
                      _slip_resistance[_qp][i];
                      
      if (_creep_activated) {

        if (_creep_resistance_function) {
          creep_resistance = _creep_resistance_function->value(_t, _q_point[_qp]);
          stress_ratio = std::abs(effective_stress / creep_resistance);
        } else {
          creep_resistance = _slip_resistance[_qp][i];
        }

        dslip_dtau[i] += creep_ao / _creep_xm * std::pow(stress_ratio, 1.0 / _creep_xm - 1.0) /
                         creep_resistance * tertiary_creep;
      }
    }
  }
}

bool
CrystalPlasticityDislocationUpdateNi::areConstitutiveStateVariablesConverged()
{
  bool rho_converged = isConstitutiveStateVariableConverged(_rho_ssd[_qp], // NEW
                                              _rho_ssd_before_update,
                                              _previous_substep_rho_ssd,
                                              _rho_tol);
  bool r_g_prime_converged = isConstitutiveStateVariableConverged(_r_eff_g_prime[_qp], // NEW
                                              _r_eff_g_prime_before_update,
                                              _previous_substep_r_eff_g_prime,
                                              _r_prep_tol); // TODO: Needs its own tolerance
  bool r_g_pp_converged = isConstitutiveStateVariableConverged(_r_eff_g_pp[_qp], // NEW
                                              _r_eff_g_pp_before_update,
                                              _previous_substep_r_eff_g_pp,
                                              _r_prep_tol); // TODO: Needs its own tolerance
  
  return rho_converged && r_g_prime_converged && r_g_pp_converged;
}

void
CrystalPlasticityDislocationUpdateNi::updateSubstepConstitutiveVariableValues()
{
  // Update temporary variables at the end of the substep
  _previous_substep_rho_ssd = _rho_ssd[_qp];
  _previous_substep_rho_gnd_edge = _rho_gnd_edge[_qp];
  _previous_substep_rho_gnd_screw = _rho_gnd_screw[_qp];
  _previous_substep_backstress_inter = _backstress_inter[_qp]; // RENAMED
  _previous_substep_backstress_intra = _backstress_intra[_qp]; // NEW
  _previous_substep_r_eff_g_prime = _r_eff_g_prime[_qp]; // NEW
  _previous_substep_r_eff_g_pp = _r_eff_g_pp[_qp];       // NEW
}

void
CrystalPlasticityDislocationUpdateNi::cacheStateVariablesBeforeUpdate()
{
  _rho_ssd_before_update = _rho_ssd[_qp];
  _rho_gnd_edge_before_update = _rho_gnd_edge[_qp];
  _rho_gnd_screw_before_update = _rho_gnd_screw[_qp];
  _backstress_inter_before_update = _backstress_inter[_qp]; // RENAMED
  _backstress_intra_before_update = _backstress_intra[_qp]; // NEW
  _r_eff_g_prime_before_update = _r_eff_g_prime[_qp]; // NEW
  _r_eff_g_pp_before_update = _r_eff_g_pp[_qp];       // NEW
}

void
CrystalPlasticityDislocationUpdateNi::calculateStateVariableEvolutionRateComponent()
{
  Real rho_sum;
  Real T = _temperature[_qp];  // NEW
  Real T_ref = _reference_temperature;  // NEW
  
  // When T increases, the y_c increases for higher annihilation.
  Real y_c_eff = _y_c * std::exp(-(_Q_drv / _R_gas_constant) * (1.0 / T - 1.0 / T_ref));  // Q_drv needs to be calibrated for different temperatures.

  // SSD dislocation density increment
  for (const auto i : make_range(_number_slip_systems))
  {
    
    rho_sum = _rho_ssd[_qp][i] + std::abs(_rho_gnd_edge[_qp][i]) + std::abs(_rho_gnd_screw[_qp][i]);

    // Multiplication and annihilation
	// note that _slip_increment here is the rate
	// and the rate equation gets multiplied by time step in updateStateVariables
    _rho_ssd_increment[i] = _k_0 * sqrt(rho_sum) - 2 * y_c_eff * _rho_ssd[_qp][i]; // NEW
    _rho_ssd_increment[i] *= std::abs(_slip_increment[_qp][i]) / _burgers_vector_mag;

   // Implement Precipitate Softening Evolution (Priority 2) ---
    // Implement decay rule: dr/dt = -C_shear * r * |slip_rate|
   // To be considered: The tempearture effect should be within the thermally-activated
   // slip rate, which require ao and xm to be temperature-dependent. (Now ao is not) 

    // Get current slip rate 

    // Calculate g' (L12) radius evolution
    Real r_current_gp = _r_eff_g_prime[_qp][i];
    _r_eff_g_prime_increment[i] = -_C_shear_g_prime * r_current_gp * std::abs(_slip_increment[_qp][i]);

    // Calculate g'' (DO22) radius evolution
    Real r_current_gpp = _r_eff_g_pp[_qp][i];
    _r_eff_g_pp_increment[i] = -_C_shear_g_pp * r_current_gpp * std::abs(_slip_increment[_qp][i]);

  }
  
  // GND dislocation density increment
  if (_include_slip_gradients) {
    for (const auto i : make_range(_number_slip_systems))
    {

    _rho_gnd_edge_increment[i] = (-1.0) * _dslip_increment_dedge[_qp](i) / _burgers_vector_mag;
    _rho_gnd_screw_increment[i] = _dslip_increment_dscrew[_qp](i) / _burgers_vector_mag;

    }
  }
  
  // backstress increment
  ArmstrongFrederickBackstressUpdate();
}

// Armstrong-Frederick update of the backstress
void
CrystalPlasticityDislocationUpdateNi::ArmstrongFrederickBackstressUpdate()
{
// --- Pre-calculate shared values ---
  Real G_T = _shear_modulus * (1.0 - _dshear_dT * (_temperature[_qp] - _reference_temperature));
  Real b = _burgers_vector_mag;
  Real f_gpp = _f_vol_g_pp;        // Precipitate volume fraction
  Real r_init_gpp = _init_r_eff_g_pp; // Initial precipitate radius

  for (const auto i : make_range(_number_slip_systems)) 
  {

    // --- Part 1: Intergranular Backstress (Existing Armstrong-Frederick Model) ---
    _backstress_inter_increment[i] = _h * _slip_increment[_qp][i];
    _backstress_inter_increment[i] -= _h_D * _backstress_inter[_qp][i] * std::abs(_slip_increment[_qp][i]);

// --- Part 2: Intragranular Backstress (Agaram et al. Eq 9, 10, 30) --- 
    
    // Calculate L_ppt (mean free path) from Agaram Eq. (30) 
    Real r_eff_gpp = _r_eff_g_pp[_qp][i]; // Current (evolving) radius ISV
    Real L_ppt = 0.0;
    if (f_gpp > 1.0e-12 && r_init_gpp > 1.0e-12)
    {
        // Note: (pi / 3f)^(1/3)
        Real geom_factor = std::pow(M_PI / (3.0 * f_gpp), 1.0/3.0);
        L_ppt = 2.0 * (geom_factor * r_init_gpp - r_eff_gpp);
    }

    // Calculate Z_o (precipitate dislocation generation) from Agaram Eq. (3) 
    Real Z_o = 0.0;
    if (L_ppt > 1.0e-12) // Avoid division by zero
    {
        Z_o = (1.0 / b) * (_k_D / L_ppt);
    }

    // Calculate total dislocation density (rho_total)
    Real rho_total = _rho_ssd[_qp][i] + std::abs(_rho_gnd_edge[_qp][i]) + std::abs(_rho_gnd_screw[_qp][i]);
    Real k1_sqrt_rho = _k_0 * std::sqrt(rho_total); // _k_0 is our k_1

    // Calculate eta (precipitate dislocation fraction) from Agaram Eq. (10) 
    Real eta = 0.0;
    Real denominator = Z_o + k1_sqrt_rho; // Simplified form
    if (std::abs(denominator) > 1.0e-9)
    {
        eta = (_k_32 * Z_o) / denominator;
    }

    // Calculate intra-backstress evolution from Agaram Eq. (9) [cite: 158-159]
    Real tau_rss = _tau[_qp][i]; // Current resolved shear stress
     Real chi_intra = _backstress_intra[_qp][i]; // Current intra-backstress
    
    // sign(tau - chi_intra)
    Real sign_term = (tau_rss - chi_intra > 0.0) ? 1.0 : -1.0; 
    
    // Term 1: [eta * b * mu * sqrt(rho) * sign(..)]
    Real term1 = eta * b * G_T * std::sqrt(rho_total) * sign_term;
    
    // Term 2: [- chi_intra]
    Real term2 = -chi_intra;

    // d(chi_intra)/dt = k_52 * (Term1 + Term2) * |slip_rate|
    _backstress_intra_increment[i] = _k_52 * (term1 + term2) * std::abs(_slip_increment[_qp][i]);
  }
}
bool
CrystalPlasticityDislocationUpdateNi::updateStateVariables()
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
  // Backstress (Intergranular): can be both positive or negative
  for (const auto i : make_range(_number_slip_systems))
  { 
    _backstress_inter_increment[i] *= _substep_dt; // RENAMED
    _backstress_inter[_qp][i] = _previous_substep_backstress_inter[i] + _backstress_inter_increment[i]; // RENAMED
  }

  // Backstress (Intragranular): can be both positive or negative
  for (const auto i : make_range(_number_slip_systems))
  { 
    _backstress_intra_increment[i] *= _substep_dt; // NEW
    _backstress_intra[_qp][i] = _previous_substep_backstress_intra[i] + _backstress_intra_increment[i]; // NEW
  }
  // Update precipitate radii ISVs 
  for (const auto i : make_range(_number_slip_systems))
  { 
    _r_eff_g_prime_increment[i] *= _substep_dt; // NEW
    _r_eff_g_prime[_qp][i] = _previous_substep_r_eff_g_prime[i] + _r_eff_g_prime_increment[i]; // NEW
    // Add safety check for radius (cannot be negative) 
    if (_r_eff_g_prime[_qp][i] < 0.0) _r_eff_g_prime[_qp][i] = 0.0; // NEW

    _r_eff_g_pp_increment[i] *= _substep_dt; // NEW
    _r_eff_g_pp[_qp][i] = _previous_substep_r_eff_g_pp[i] + _r_eff_g_pp_increment[i]; // NEW
    if (_r_eff_g_pp[_qp][i] < 0.0) _r_eff_g_pp[_qp][i] = 0.0; // NEW
  }  
  return true;
}
