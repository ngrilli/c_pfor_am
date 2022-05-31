// Nicolò Grilli
// University of Bristol
// 12 Maggio 2022

// This is a dislocation-based constitutive model
// for ferritic reactor pressure vessel steel
// as published in:
// Ghiath Monnet, Ludovic Vincent, Lionel Gelebart
// Multiscale modeling of crystal plasticity in Reactor Pressure Vessel
// steels: Prediction of irradiation hardening
// Journal of Nuclear Materials 514 (2019) 128-138
// Additionally GND density due to slip gradients
// are also included
// This model is meant to be used with BCC slip systems only

#include "CrystalPlasticityIrradiatedRPVSteel.h"
#include "libmesh/int_range.h"
#include <cmath>

registerMooseObject("TensorMechanicsApp", CrystalPlasticityIrradiatedRPVSteel);

InputParameters
CrystalPlasticityIrradiatedRPVSteel::validParams()
{
  InputParameters params = CrystalPlasticityDislocationUpdateBase::validParams();
  params.addClassDescription("Dislocation based model for crystal plasticity "
                             "using the stress update code."
							 "Irradiation damage in RPV steel.");
							 
  params.addParam<Real>("burgers_vector_mag",0.000256,"Magnitude of the Burgers vector (micron)");
  params.addParam<Real>("shear_modulus",86000.0,"Shear modulus in Taylor hardening law G");
  params.addParam<Real>("RT_shear_modulus",86000.0,"Shear modulus at room temperature");
  params.addParam<Real>("a_self",0.1,"Self interaction coefficient of the slip systems");
  params.addParam<Real>("a_col",0.7,"Collinear interaction coefficient of the slip systems");
  params.addParam<Real>("K_Hall_Petch",480.0,"Hall-Petch effect prefactor (MPa micron^(1/2))");
  params.addParam<Real>("d_grain",6.9,"Average grain size (micron)");
  params.addParam<Real>("rho_carbide",0.0608,"Carbide planar density (micron)^{-2}");
  params.addParam<Real>("a_carbide",0.0,"Carbide interaction coefficient");
  params.addParam<Real>("C_DL_diameter",0.0256,"Average diameter of irradiation dislocation loops (micron)");
  params.addParam<Real>("a_DL",0.0,"Irradiation dislocation loops interaction coefficient");
  params.addParam<Real>("C_SC_diameter",0.0256,"Average diameter of irradiation solute clusters (micron)");
  params.addParam<Real>("a_SC",0.0,"Solute clusters interaction coefficient");
  params.addParam<Real>("rho_ref",1.0,"Reference dislocation density at which the interaction "
                                      "matrix between slip system is the reference matrix "
									  "(1 / micron^2)");
  
  
  

							 
  // TO DO - add or remove parameters						 
  params.addParam<Real>("ao", 0.001, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");  

  params.addParam<Real>("alpha_0",0.3,"Prefactor of Taylor hardening law, alpha");
  params.addParam<Real>("r", 1.4, "Latent hardening coefficient");
  params.addParam<Real>("tau_c_0", 0.112, "Peierls stress");
  params.addParam<Real>("k_0",100.0,"Coefficient K in SSD evolution, representing accumulation rate");
  params.addParam<Real>("y_c",0.0026,"Critical annihilation diameter");
  params.addParam<Real>("rho_tol",1.0,"Tolerance on dislocation density update");
  
  
  
  params.addParam<Real>("init_rho_ssd",10.0,"Initial dislocation density (micron)^{-2}");
  params.addParam<Real>("init_rho_gnd_edge",0.0,"Initial dislocation density (micron)^{-2}");
  params.addParam<Real>("init_rho_gnd_screw",0.0,"Initial dislocation density (micron)^{-2}");
  params.addParam<Real>("init_C_DL",0.0,"Initial concentration of irradiation dislocation loops");
  params.addParam<Real>("init_C_SC",0.0,"Initial concentration of irradiation solute clusters");
  params.addParam<MaterialPropertyName>(
      "total_twin_volume_fraction",
      "Total twin volume fraction, if twinning is considered in the simulation");
  params.addParam<UserObjectName>("read_initial_gnd_density",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element value "
                                  "of the initial GND density");
  params.addCoupledVar("dslip_increment_dedge",0.0,"Directional derivative of the slip rate along the edge motion direction.");
  params.addCoupledVar("dslip_increment_dscrew",0.0,"Directional derivative of the slip rate along the screw motion direction.");
  return params;
}

CrystalPlasticityIrradiatedRPVSteel::CrystalPlasticityIrradiatedRPVSteel(
    const InputParameters & parameters)
  : CrystalPlasticityDislocationUpdateBase(parameters),
  
    // Constitutive model parameters
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")),
	_shear_modulus(getParam<Real>("shear_modulus")),
	_RT_shear_modulus(getParam<Real>("RT_shear_modulus")),
    _a_self(getParam<Real>("a_self")),
	_a_col(getParam<Real>("a_col")),
    _K_Hall_Petch(getParam<Real>("K_Hall_Petch")),
	_d_grain(getParam<Real>("d_grain")),
	_rho_carbide(getParam<Real>("rho_carbide")),
	_a_carbide(getParam<Real>("a_carbide")),
	_C_DL_diameter(getParam<Real>("C_DL_diameter")),
    _a_DL(getParam<Real>("a_DL")),
	_C_SC_diameter(getParam<Real>("C_SC_diameter")),
    _a_SC(getParam<Real>("a_SC")),
	_rho_ref(getParam<Real>("rho_ref")),
	
	
	
	
	
	
	
	// TO DO
    _ao(getParam<Real>("ao")),
    _xm(getParam<Real>("xm")),

	_alpha_0(getParam<Real>("alpha_0")),
    _r(getParam<Real>("r")),
	_tau_c_0(getParam<Real>("tau_c_0")),
	_k_0(getParam<Real>("k_0")),
	_y_c(getParam<Real>("y_c")),
	
	
	
	
	
	
	// Initial values of the state variables
    _init_rho_ssd(getParam<Real>("init_rho_ssd")),
    _init_rho_gnd_edge(getParam<Real>("init_rho_gnd_edge")),
    _init_rho_gnd_screw(getParam<Real>("init_rho_gnd_screw")),
	
	// Initial values of the irradiation defects
	_init_C_DL(getParam<Real>("init_C_DL")),
	_init_C_SC(getParam<Real>("init_C_SC")),
	
	// Tolerance on dislocation density update
	_rho_tol(getParam<Real>("rho_tol")),
	
    // state variables
	
	// _rho_ssd corresponds to rho_s in equation (18)
    _rho_ssd(declareProperty<std::vector<Real>>("rho_ssd")),
    _rho_ssd_old(getMaterialPropertyOld<std::vector<Real>>("rho_ssd")),
	
	// GND dislocation densities: not in the original model
    _rho_gnd_edge(declareProperty<std::vector<Real>>("rho_gnd_edge")),
   	_rho_gnd_edge_old(getMaterialPropertyOld<std::vector<Real>>("rho_gnd_edge")),
  	_rho_gnd_screw(declareProperty<std::vector<Real>>("rho_gnd_screw")),
    _rho_gnd_screw_old(getMaterialPropertyOld<std::vector<Real>>("rho_gnd_screw")),

    // C_DL: concentration of dislocation loops induced by irradiation
	// on each slip system
	_C_DL(declareProperty<std::vector<Real>>("C_DL")),
	_C_DL_old(getMaterialPropertyOld<std::vector<Real>>("C_DL")),
	
	// C_SC: concentration of solute clusters
	// on each slip system
	_C_SC(declareProperty<std::vector<Real>>("C_SC")),
	_C_SC_old(getMaterialPropertyOld<std::vector<Real>>("C_SC")),
	
	// increment of state variables
    _rho_ssd_increment(_number_slip_systems, 0.0),
    _rho_gnd_edge_increment(_number_slip_systems, 0.0),
    _rho_gnd_screw_increment(_number_slip_systems, 0.0),
	_C_DL_increment(_number_slip_systems, 0.0),
	_C_SC_increment(_number_slip_systems, 0.0),
	
	// resize local caching vectors used for substepping
    _previous_substep_rho_ssd(_number_slip_systems, 0.0),
	_previous_substep_rho_gnd_edge(_number_slip_systems, 0.0),
	_previous_substep_rho_gnd_screw(_number_slip_systems, 0.0),
	_previous_substep_C_DL(_number_slip_systems, 0.0),
	_previous_substep_C_SC(_number_slip_systems, 0.0),
	
    _rho_ssd_before_update(_number_slip_systems, 0.0),
    _rho_gnd_edge_before_update(_number_slip_systems, 0.0),
    _rho_gnd_screw_before_update(_number_slip_systems, 0.0),
    _C_DL_before_update(_number_slip_systems, 0.0),
	_C_SC_before_update(_number_slip_systems, 0.0),

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

    // store edge and screw slip directions to calculate directional derivatives
    // of the plastic slip rate	
    _edge_slip_direction(declareProperty<std::vector<Real>>("edge_slip_direction")),
	_screw_slip_direction(declareProperty<std::vector<Real>>("screw_slip_direction")),
	
	// Total density of dislocations including SSD and GND on each slip system
	_rho_tot(_number_slip_systems, 0.0),
	
	// Total density of local obstacles for each slip system
	_rho_obstacles(_number_slip_systems, 0.0),
	
	// Average obstacles strength for each slip system
	// according to equation (6)
	_obstacles_strength(_number_slip_systems, 0.0),
	
	// Self interaction stress tau_self for each slip system
	_tau_self(_number_slip_systems, 0.0),
	
	// Line tension slip resistance for each slip system
	_tau_line_tension(_number_slip_systems, 0.0),
	
	// Reference interaction matrix between slip systems
	_a_ref(_number_slip_systems, _number_slip_systems),
	
	// Corrected interaction matrix between slip systems
	// that accounts for the logarithmic correction in equation (7)
	_a_slip_slip_interaction(_number_slip_systems, _number_slip_systems)
{
}

void
CrystalPlasticityIrradiatedRPVSteel::initQpStatefulProperties()
{
  // Slip resistance is resized here
  // _slip_resistance is tau_c in equation (11)
  CrystalPlasticityDislocationUpdateBase::initQpStatefulProperties();
  
  Real taylor_hardening;

  // Initialize and dislocation density size
  _rho_ssd[_qp].resize(_number_slip_systems);
  _rho_gnd_edge[_qp].resize(_number_slip_systems);
  _rho_gnd_screw[_qp].resize(_number_slip_systems);
  
  // Initialize irradiation defects size
  _C_DL[_qp].resize(_number_slip_systems);
  _C_SC[_qp].resize(_number_slip_systems);
  
  // Initialize dislocation densities
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
  }
  
  // Initialize irradiation defects
  for (const auto i : make_range(_number_slip_systems)) 
  {
    _C_DL[_qp][i] = _init_C_DL;
	_C_SC[_qp][i] = _init_C_SC;
  }

  // Initialize constant reference interaction matrix between slip systems
  initializeReferenceInteractionMatrix();
  
  // Function to calculate slip resistance can be called
  // here because the state variables are already assigned
  calculateSlipResistance();
  

  
  // TO DO
  
  // Initialize value of the slip resistance
  // as a function of the dislocation density
  // and irradiation defect density
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
	                          * std::sqrt(taylor_hardening));
	
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

// Initialize constant reference interaction matrix between slip systems
// based on Figure 1
// This is meant to work for BCC crystals only
void
CrystalPlasticityIrradiatedRPVSteel::initializeReferenceInteractionMatrix()
{
  // Initially fill with self interaction
  for (const auto i : make_range(_number_slip_systems)) {
    for (const auto j : make_range(_number_slip_systems)) {
		
      _a_ref(i,j) = _a_self;	
		
	}  
  }

  // Fill the collinear in the 3x3 diagonal blocks
  for (unsigned int i = 0; i < (_number_slip_systems/3); ++i) {
	  
    _a_ref(3*i,3*i+1) = _a_col;
	_a_ref(3*i+1,3*i) = _a_col;
    _a_ref(3*i,3*i+2) = _a_col;
	_a_ref(3*i+2,3*i) = _a_col;
    _a_ref(3*i+1,3*i+2) = _a_col;
    _a_ref(3*i+2,3*i+1) = _a_col;	

  }
	
  // Fill the collinear in the 3x3 cross slip blocks
  // Assigning one by one because there is little regularity
  _a_ref(0,13) = _a_col;
  _a_ref(1,14) = _a_col;
  _a_ref(2,12) = _a_col;
  _a_ref(3,16) = _a_col;
  _a_ref(4,15) = _a_col;
  _a_ref(5,17) = _a_col;
  _a_ref(6,20) = _a_col;
  _a_ref(7,18) = _a_col;
  _a_ref(8,19) = _a_col;
  _a_ref(9,22) = _a_col;
  _a_ref(10,23) = _a_col;
  _a_ref(11,21) = _a_col;
  
  // add the transpose part
  _a_ref(13,0) = _a_col;
  _a_ref(14,1) = _a_col;
  _a_ref(12,2) = _a_col;
  _a_ref(16,3) = _a_col;
  _a_ref(15,4) = _a_col;
  _a_ref(17,5) = _a_col;
  _a_ref(20,6) = _a_col;
  _a_ref(18,7) = _a_col;
  _a_ref(19,8) = _a_col;
  _a_ref(22,9) = _a_col;
  _a_ref(23,10) = _a_col;
  _a_ref(21,11) = _a_col;
  
}

// Logarithmic correction to the interaction matrix in equation (7)
void
CrystalPlasticityIrradiatedRPVSteel::logarithmicCorrectionInteractionMatrix()
{
  // rho obstacles must be already calculated at this point
  
  // temporary variable to store the logarithmic correction prefactor
  Real temp_log_factor;
  
  for (const auto i : make_range(_number_slip_systems)) {
    for (const auto j : make_range(_number_slip_systems)) {
	
      temp_log_factor = std::log(0.35 * _burgers_vector_mag * std::sqrt(_rho_obstacles[i]));
	  temp_log_factor /= std::log(0.35 * _burgers_vector_mag * std::sqrt(_rho_ref));
	  temp_log_factor *= 0.8;
	  temp_log_factor += 0.2;
  
      _a_slip_slip_interaction(i,j) = temp_log_factor * temp_log_factor * _a_ref(i,j);
  
	}
  }
}

// Calculate Schmid tensor and
// store edge and screw slip directions to calculate directional derivatives
// of the plastic slip rate
// this is a general function and not specific to the constitutive model
void
CrystalPlasticityIrradiatedRPVSteel::calculateSchmidTensor(
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
CrystalPlasticityIrradiatedRPVSteel::setInitialConstitutiveVariableValues()
{
  // Would also set old dislocation densities here if included in this model
  _rho_ssd[_qp] = _rho_ssd_old[_qp];
  _previous_substep_rho_ssd = _rho_ssd_old[_qp];
  _rho_gnd_edge[_qp] = _rho_gnd_edge_old[_qp];
  _previous_substep_rho_gnd_edge = _rho_gnd_edge_old[_qp];
  _rho_gnd_screw[_qp] = _rho_gnd_screw_old[_qp];
  _previous_substep_rho_gnd_screw = _rho_gnd_screw_old[_qp];
  _C_DL[_qp] = _C_DL_old[_qp];
  _previous_substep_C_DL = _C_DL_old[_qp];
  _C_SC[_qp] = _C_SC_old[_qp];
  _previous_substep_C_SC = _C_SC_old[_qp];
}

void
CrystalPlasticityIrradiatedRPVSteel::setSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _rho_ssd[_qp] = _previous_substep_rho_ssd;
  _rho_gnd_edge[_qp] = _previous_substep_rho_gnd_edge;
  _rho_gnd_screw[_qp] = _previous_substep_rho_gnd_screw;
  _C_DL[_qp] = _previous_substep_C_DL;
  _C_SC[_qp] = _previous_substep_C_SC;
}

// Slip resistance can be calculated from dislocation density here only
// because it is the first method in which it is used,
// while calculateConstitutiveSlipDerivative is called afterwards
bool
CrystalPlasticityIrradiatedRPVSteel::calculateSlipRate()
{
  calculateSlipResistance();
	
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_increment[_qp][i] =
        _ao * std::pow(std::abs(_tau[_qp][i] / _slip_resistance[_qp][i]), 1.0 / _xm);
    if (_tau[_qp][i] < 0.0)
      _slip_increment[_qp][i] *= -1.0;

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

// Slip resistance based on equation (11)
void
CrystalPlasticityIrradiatedRPVSteel::calculateSlipResistance()
{
  // Calculate total density of local obstacles
  calculateObstaclesDensity();
  
  // Calculate logarithmic correction to slip-slip
  // interaction matrix
  logarithmicCorrectionInteractionMatrix();
  
  // Calculate average obstacles strength
  calculateObstaclesStrength();
	
  // Calculate different contributions to the CRSS
  calculateSelfInteractionSlipResistance();
  calculateHallPetchSlipResistance();
  calculateLineTensionSlipResistance();

  // sum the contributions to the CRSS
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance[_qp][i] = _tau_self[i];
	_slip_resistance[_qp][i] += _tau_Hall_Petch;
	_slip_resistance[_qp][i] += _tau_line_tension[i];
  }  
  

  // TO DO
  Real taylor_hardening;
	
	
  // TO DO	
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
	                          * std::sqrt(taylor_hardening));
	
  }

}

// Calculate total density of local obstacles based on equation (5)
// this function stores at the same time the sum of SSD and GND
// which is needed later in the obstacles strength calculation
void
CrystalPlasticityIrradiatedRPVSteel::calculateObstaclesDensity()
{
  // Define total dislocation density
  // including GND and SSD on each slip system 
	
  for (const auto i : make_range(_number_slip_systems))	
  {
    _rho_tot[i] = _rho_ssd[_qp][i] 
                + std::abs(_rho_gnd_edge[_qp][i])
                + std::abs(_rho_gnd_screw[_qp][i]);
  }
	
  for (const auto i : make_range(_number_slip_systems)) 
  {
    _rho_obstacles[i] = 0.0;
		
    for (const auto j : make_range(_number_slip_systems)) 
    {
		 
      if (j != i) {

	    _rho_obstacles[i] += _rho_tot[j];
		
	  } 
    }
  }		
  
  // Add carbide density
  for (const auto i : make_range(_number_slip_systems)) {
	  
    _rho_obstacles[i] += _rho_carbide;
	  
  }
  
  // Add irradiation defects densities
  for (const auto i : make_range(_number_slip_systems)) {
	  
    _rho_obstacles[i] += _C_DL_diameter * _C_DL[_qp][i];
	_rho_obstacles[i] += _C_SC_diameter * _C_SC[_qp][i];
	  
  }
  
}

// Calculate average obstacle strength based on equation (6)
void
CrystalPlasticityIrradiatedRPVSteel::calculateObstaclesStrength()
{
  // temporary variable to sum the contributions
  // of different mechanisms for different slip systems
  Real temp_sqrt_argument;
  
  for (const auto i : make_range(_number_slip_systems)) {
	  
    temp_sqrt_argument = 0.0;
	  
    // Other slip systems contribution 
    for (const auto j : make_range(_number_slip_systems)) {
		
      if (j != i) {

	    temp_sqrt_argument += _a_slip_slip_interaction(i,j) * _rho_tot[j];		  
		  
	  }	
	}
	
	// carbide contribution
	temp_sqrt_argument += _a_carbide * _rho_carbide;

    // irradiation damage contribution
    // would be better to store rho_DL and rho_SC
    // into state variables because the multiplication is repeated
    temp_sqrt_argument += _a_DL * _C_DL_diameter * _C_DL[_qp][i];
    temp_sqrt_argument += _a_SC	* _C_SC_diameter * _C_SC[_qp][i];

    _obstacles_strength[i] = std::sqrt(temp_sqrt_argument) / std::sqrt(_rho_obstacles[i]);
	
  }

}

// Self-interaction slip resistance based on equation (12)
void
CrystalPlasticityIrradiatedRPVSteel::calculateSelfInteractionSlipResistance()
{
  for (const auto i : make_range(_number_slip_systems))	{
	  
    _tau_self[i] = _burgers_vector_mag * _shear_modulus 
	             * std::sqrt(_a_self * _rho_tot[i]);
	
  }

}

// Hall-Petch slip resistance based on equation (13)
// note that if GND are activated, part of the Hall-Petch effect
// will be provided by the GNDs
void
CrystalPlasticityIrradiatedRPVSteel::calculateHallPetchSlipResistance()
{
	
  _tau_Hall_Petch = (_shear_modulus / _RT_shear_modulus) 
                  * (_K_Hall_Petch / std::sqrt(_d_grain));

}

// Line tension slip resistance based on equation (15)
// note that if GND are activated, part of the Hall-Petch effect
// will be provided by the GNDs
void
CrystalPlasticityIrradiatedRPVSteel::calculateLineTensionSlipResistance()
{
  for (const auto i : make_range(_number_slip_systems))	{
	  
	// TO DO: add density and alpha^s
    _tau_line_tension[i] = _shear_modulus * _burgers_vector_mag;
	
  }
  
}

void
CrystalPlasticityIrradiatedRPVSteel::calculateEquivalentSlipIncrement(
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
CrystalPlasticityIrradiatedRPVSteel::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  for (const auto i : make_range(_number_slip_systems))
  {
    if (MooseUtils::absoluteFuzzyEqual(_tau[_qp][i], 0.0))
      dslip_dtau[i] = 0.0;
    else
      dslip_dtau[i] = _ao / _xm *
                      std::pow(std::abs(_tau[_qp][i] / _slip_resistance[_qp][i]), 1.0 / _xm - 1.0) /
                      _slip_resistance[_qp][i];
  }
}

bool
CrystalPlasticityIrradiatedRPVSteel::areConstitutiveStateVariablesConverged()
{
  return isConstitutiveStateVariableConverged(_rho_ssd[_qp],
                                              _rho_ssd_before_update,
                                              _previous_substep_rho_ssd,
                                              _rho_tol);

  // How do we check the tolerance of GNDs and is it needed?
											  
}

void
CrystalPlasticityIrradiatedRPVSteel::updateSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _previous_substep_rho_ssd = _rho_ssd[_qp];
  _previous_substep_rho_gnd_edge = _rho_gnd_edge[_qp];
  _previous_substep_rho_gnd_screw = _rho_gnd_screw[_qp];
  _previous_substep_C_DL = _C_DL[_qp];
  _previous_substep_C_SC = _C_SC[_qp];
}

void
CrystalPlasticityIrradiatedRPVSteel::cacheStateVariablesBeforeUpdate()
{
  _rho_ssd_before_update = _rho_ssd[_qp];
  _rho_gnd_edge_before_update = _rho_gnd_edge[_qp];
  _rho_gnd_screw_before_update = _rho_gnd_screw[_qp];
  _C_DL_before_update = _C_DL[_qp];
  _C_SC_before_update = _C_SC[_qp];
}

void
CrystalPlasticityIrradiatedRPVSteel::calculateStateVariableEvolutionRateComponent()
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
}

bool
CrystalPlasticityIrradiatedRPVSteel::updateStateVariables()
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
  
  // GND edge: can be both positive and negative
  for (const auto i : make_range(_number_slip_systems))
  { 
    _rho_gnd_edge_increment[i] *= _substep_dt;
    _rho_gnd_edge[_qp][i] = _previous_substep_rho_gnd_edge[i] + _rho_gnd_edge_increment[i];
  }
  
  // GND screw: can be both positive and negative
  for (const auto i : make_range(_number_slip_systems))
  { 
    _rho_gnd_screw_increment[i] *= _substep_dt;
    _rho_gnd_screw[_qp][i] = _previous_substep_rho_gnd_screw[i] + _rho_gnd_screw_increment[i];
  }
  
  return true;
}
