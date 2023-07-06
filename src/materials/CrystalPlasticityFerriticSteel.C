// Michael Salvini
// Nicolò Grilli
// University of Bristol
// 8 Maggio 2023

// This is a dislocation-based constitutive model
// for ferritic reactor pressure vessel steel.
// Additionally GND density due to slip gradients
// are also included
// This model is meant to be used with BCC slip systems only

#include "CrystalPlasticityFerriticSteel.h"
#include "libmesh/int_range.h"
#include <cmath>

registerMooseObject("TensorMechanicsApp", CrystalPlasticityFerriticSteel);

InputParameters
CrystalPlasticityFerriticSteel::validParams()
{
  InputParameters params = CrystalPlasticityDislocationUpdateBase::validParams();
  params.addClassDescription("Dislocation based model for crystal plasticity "
                             "using the stress update code. "
							 "Irradiation hardening in RPV steel with slip gradients. "
							 "This model is meant to be used with BCC slip systems only. ");
  params.addParam<Real>("burgers_vector_mag",0.000248,"Magnitude of the Burgers vector (micron)");
  params.addParam<Real>("shear_modulus",86000.0,"Shear modulus in Taylor hardening law G");
  params.addParam<Real>("a_self",0.1,"Self interaction coefficient of the slip systems");
  params.addParam<Real>("a_col",0.7,"Collinear interaction coefficient of the slip systems");
  
  // Irradiation defects parameters
  params.addParam<Real>("C_DL_diameter",0.0256,"Average diameter of irradiation dislocation loops (micron)");
  params.addParam<Real>("a_DL",0.25,"Irradiation dislocation loops interaction coefficient");
  params.addParam<Real>("C_SC_diameter",0.0256,"Average diameter of irradiation solute clusters (micron)");
  params.addParam<Real>("a_SC",0.04,"Solute clusters interaction coefficient");
  params.addParam<Real>("lambda_DL", 1.0,"prefactor of the irradiation dislocation loops evolution law (adimensional)");
  params.addParam<Real>("lambda_SC", 1.0,"prefactor of the irradiation solute cluster evolution law (adimensional)");
  params.addParam<bool>("is_irradiated", false, "Skip irradiation calculations if false.");
  
  // Slip and creep rate parameters
  params.addParam<Real>("ao", 0.00001, "slip rate coefficient (s^{-1})");
  params.addParam<Real>("xm", 0.01, "exponent for slip rate");
  params.addParam<bool>("creep_activated", false, "Activate creep strain rate.");
  params.addParam<Real>("creep_ao", 0.0, "creep rate coefficient");
  params.addParam<Real>("creep_xm", 0.1, "exponent for creep rate");
  params.addParam<Real>("m_exponent", 0.0, "Exponent on time in power-law equation");
  params.addParam<Real>("max_stress_ratio", 10.0, "Maximum ratio between RSS and CRSS above which slip law becomes linear");
  params.addParam<Real>("reduced_ao", 0.00001, "slip rate coefficient (s^{-1}) of the linear law once max stress ratio is exceeded");

  // Constant slip resistances of the slip systems
  // Hall-Petch effect must be included in these constants
  // Carbide-dislocation interaction must be included in these constants
  params.addParam<Real>("const_slip_resistance_110", 360.0, "Constant slip resistance of 110 slip planes (MPa)");
  params.addParam<Real>("const_slip_resistance_112_TW", 410.0, "Constant slip resistance of 112 slip planes in twinning direction (MPa)");
  params.addParam<Real>("const_slip_resistance_112_AT", 480.0, "Constant slip resistance of 112 slip planes in anti-twinning direction (MPa)");
  
  // Dislocation multiplication and annihilation parameters
  params.addParam<Real>("k_0",1.0,"Coefficient in SSD evolution, representing the constant accumulation rate");  
  params.addParam<Real>("y_c",0.0026,"Critical annihilation diameter");
  
  // Initial values of the state variables
  params.addParam<Real>("init_rho_ssd",10.0,"Initial dislocation density (micron)^{-2}");
  params.addParam<Real>("init_rho_gnd_edge",0.0,"Initial dislocation density (micron)^{-2}");
  params.addParam<Real>("init_rho_gnd_screw",0.0,"Initial dislocation density (micron)^{-2}");
  params.addParam<Real>("init_C_DL",0.0,"Initial concentration of irradiation dislocation loops");
  params.addParam<Real>("init_C_SC",0.0,"Initial concentration of irradiation solute clusters");
  params.addParam<Real>("rho_tol",1.0,"Tolerance on dislocation density update");
  
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

CrystalPlasticityFerriticSteel::CrystalPlasticityFerriticSteel(
    const InputParameters & parameters)
  : CrystalPlasticityDislocationUpdateBase(parameters),

    // Constitutive model parameters
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")),
	_shear_modulus(getParam<Real>("shear_modulus")),
    _a_self(getParam<Real>("a_self")),
    _a_col(getParam<Real>("a_col")),
	_C_DL_diameter(getParam<Real>("C_DL_diameter")),
    _a_DL(getParam<Real>("a_DL")),
	_C_SC_diameter(getParam<Real>("C_SC_diameter")),
    _a_SC(getParam<Real>("a_SC")),
    _lambda_DL(getParam<Real>("lambda_DL")),
    _lambda_SC(getParam<Real>("lambda_SC")),
    _is_irradiated(getParam<bool>("is_irradiated")),
	_ao(getParam<Real>("ao")),
	_xm(getParam<Real>("xm")),
	_creep_activated(getParam<bool>("creep_activated")),
    _creep_ao(getParam<Real>("creep_ao")),
    _creep_xm(getParam<Real>("creep_xm")),
    _m_exponent(getParam<Real>("m_exponent")),
    _max_stress_ratio(getParam<Real>("max_stress_ratio")),
    _reduced_ao(getParam<Real>("reduced_ao")),
    _const_slip_resistance_110(getParam<Real>("const_slip_resistance_110")),
    _const_slip_resistance_112_TW(getParam<Real>("const_slip_resistance_112_TW")),
    _const_slip_resistance_112_AT(getParam<Real>("const_slip_resistance_112_AT")),
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

    // Constant slip resistance \tau_0
    // for the different slip systems
    // Hall-Petch effect must be included in these constants
    // Carbide-dislocation interaction must be included in these constants
    _const_slip_resistance(_number_slip_systems, 0.0),

	// Reference interaction matrix between slip systems
	_a_ref(_number_slip_systems, _number_slip_systems)

{
}

void
CrystalPlasticityFerriticSteel::initQpStatefulProperties()
{
  // _slip_resistance and _slip_increment are resized here
  CrystalPlasticityDislocationUpdateBase::initQpStatefulProperties();

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

  // Initialize the constant slip resistance \tau_0
  initializeConstSlipResistance();

  // Function to calculate slip resistance can be called
  // here because the state variables are already assigned
  // and _slip_resistance[_qp][i] is initialized
  calculateSlipResistance();

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
// based on Figure 1 in Monnet 2019
// This is meant to work for BCC crystals only
void
CrystalPlasticityFerriticSteel::initializeReferenceInteractionMatrix()
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

// Initialize the constant slip resistance \tau_0
// for the order of the slip systems see Table 1 in:
// BCC single crystal plasticity modeling and its
// experimental identification
// T Yalcinkaya et al 2008 Modelling Simul. Mater. Sci. Eng. 16 085007
// https://iopscience.iop.org/article/10.1088/0965-0393/16/8/085007/pdf
// Hall-Petch effect must be included in these constants
// Carbide-dislocation interaction must be included in these constants
void
CrystalPlasticityFerriticSteel::initializeConstSlipResistance() {

  // 110 slip planes
  for (const auto i : make_range(_number_slip_systems/2)) {

    _const_slip_resistance[i] = _const_slip_resistance_110;

  }

  _const_slip_resistance[12] = _const_slip_resistance_112_AT;
  _const_slip_resistance[13] = _const_slip_resistance_112_TW;
  _const_slip_resistance[14] = _const_slip_resistance_112_AT;
  _const_slip_resistance[15] = _const_slip_resistance_112_AT;
  _const_slip_resistance[16] = _const_slip_resistance_112_TW;
  _const_slip_resistance[17] = _const_slip_resistance_112_AT;
  _const_slip_resistance[18] = _const_slip_resistance_112_AT;
  _const_slip_resistance[19] = _const_slip_resistance_112_AT;
  _const_slip_resistance[20] = _const_slip_resistance_112_TW;
  _const_slip_resistance[21] = _const_slip_resistance_112_AT;
  _const_slip_resistance[22] = _const_slip_resistance_112_TW;
  _const_slip_resistance[23] = _const_slip_resistance_112_AT;
}

// Calculate Schmid tensor and
// store edge and screw slip directions to calculate directional derivatives
// of the plastic slip rate
// this is a general function and not specific to the constitutive model
void
CrystalPlasticityFerriticSteel::calculateSchmidTensor(
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
CrystalPlasticityFerriticSteel::setInitialConstitutiveVariableValues()
{
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
CrystalPlasticityFerriticSteel::setSubstepConstitutiveVariableValues()
{
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
CrystalPlasticityFerriticSteel::calculateSlipRate()
{
  // Ratio between RSS and CRSS
  // temporary variable for each slip system
  Real stress_ratio;

  // calculate and store _slip_resistance[_qp][i]
  calculateSlipResistance();
  
  for (const auto i : make_range(_number_slip_systems))
  { 
    stress_ratio = std::abs(_tau[_qp][i] / _slip_resistance[_qp][i]);
    	
    _slip_increment[_qp][i] = _ao * std::pow(stress_ratio, 1.0 / _xm);
    
    if (_creep_activated) { // add creep rate
	  _slip_increment[_qp][i] += _creep_ao * std::pow(stress_ratio, 1.0 / _creep_xm) * std::pow(_t, _m_exponent);
	}
	
	if (stress_ratio > _max_stress_ratio) { // trigger linear slip law above stress ratio threshold
      _slip_increment[_qp][i] = _ao * std::pow(_max_stress_ratio, 1.0 / _xm);
      _slip_increment[_qp][i] += _reduced_ao * (stress_ratio - _max_stress_ratio);
	}
      
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

// Slip resistance
void
CrystalPlasticityFerriticSteel::calculateSlipResistance()
{
  // Calculate total density of local obstacles
  calculateObstaclesDensity();

  // sum the contributions to the CRSS based on equation (11)
  for (const auto i : make_range(_number_slip_systems))
  {

    _slip_resistance[_qp][i] = _const_slip_resistance[i];

    _slip_resistance[_qp][i] += _burgers_vector_mag * _shear_modulus * std::sqrt(_rho_obstacles[i]);
    
  }
}

// Calculate local obstacles density _rho_obstacles
// which includes also the interaction strength
// this function stores at the same time the sum of SSD and GND: _rho_tot
void
CrystalPlasticityFerriticSteel::calculateObstaclesDensity()
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

      _rho_obstacles[i] += _a_ref(i,j) * _rho_tot[j];
	      
    }
  }

  // Add irradiation defects densities
  // and their interaction coefficients
  
  if (_is_irradiated) {
	  
    for (const auto i : make_range(_number_slip_systems)) {

      _rho_obstacles[i] += _a_DL * _C_DL_diameter * _C_DL[_qp][i];
	  _rho_obstacles[i] += _a_SC * _C_SC_diameter * _C_SC[_qp][i];

    }	    
  }

}

void
CrystalPlasticityFerriticSteel::calculateEquivalentSlipIncrement(
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
CrystalPlasticityFerriticSteel::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  // Ratio between effective stress and CRSS
  // temporary variable for each slip system
  Real stress_ratio;
	
  for (const auto i : make_range(_number_slip_systems))
  {  
	  
    if (MooseUtils::absoluteFuzzyEqual(_tau[_qp][i], 0.0)) {
		
      dslip_dtau[i] = 0.0;
      		
	} else {
		
	  stress_ratio = std::abs(_tau[_qp][i] / _slip_resistance[_qp][i]);

      dslip_dtau[i] = _ao / _xm *
                      std::pow(stress_ratio, 1.0 / _xm - 1.0) /
                      _slip_resistance[_qp][i];
                      
      if (_creep_activated) { // add creep rate
	    dslip_dtau[i] += _creep_ao / _creep_xm * std::pow(stress_ratio, 1.0 / _creep_xm - 1.0) /
	                     _slip_resistance[_qp][i] * std::pow(_t, _m_exponent);
	  }
	  
	  if (stress_ratio > _max_stress_ratio) { // trigger linear slip law above stress ratio threshold
	    dslip_dtau[i] = _reduced_ao;
	  }              
	}
  }
}

bool
CrystalPlasticityFerriticSteel::areConstitutiveStateVariablesConverged()
{
  return isConstitutiveStateVariableConverged(_rho_ssd[_qp],
                                              _rho_ssd_before_update,
                                              _previous_substep_rho_ssd,
                                              _rho_tol);
}

void
CrystalPlasticityFerriticSteel::updateSubstepConstitutiveVariableValues()
{
  _previous_substep_rho_ssd = _rho_ssd[_qp];
  _previous_substep_rho_gnd_edge = _rho_gnd_edge[_qp];
  _previous_substep_rho_gnd_screw = _rho_gnd_screw[_qp];
  _previous_substep_C_DL = _C_DL[_qp];
  _previous_substep_C_SC = _C_SC[_qp];
}

void
CrystalPlasticityFerriticSteel::cacheStateVariablesBeforeUpdate()
{
  _rho_ssd_before_update = _rho_ssd[_qp];
  _rho_gnd_edge_before_update = _rho_gnd_edge[_qp];
  _rho_gnd_screw_before_update = _rho_gnd_screw[_qp];
  _C_DL_before_update = _C_DL[_qp];
  _C_SC_before_update = _C_SC[_qp];
}

// Note that calculateStateVariableEvolutionRateComponent is called in solveStateVariables
// therefore after calculateSlipResistance
void
CrystalPlasticityFerriticSteel::calculateStateVariableEvolutionRateComponent()
{
  // Calculate increment of SSD
  calculateSSDincrement();

  // Calculate increment of irradiation defects
  if (_is_irradiated) {

    calculateDLincrement();
    calculateSCincrement();  
      
  }

  // GND dislocation density increment
  for (const auto i : make_range(_number_slip_systems))
  {

    _rho_gnd_edge_increment[i] = (-1.0) * _dslip_increment_dedge[_qp](i) / _burgers_vector_mag;
    _rho_gnd_screw_increment[i] = _dslip_increment_dscrew[_qp](i) / _burgers_vector_mag;

  }
}

// Calculate the SSD increment based on equation (18)
void
CrystalPlasticityFerriticSteel::calculateSSDincrement()
{
  // Multiplication and annihilation
  // note that _slip_increment here is the rate
  // and _rho_ssd_increment is also the rate
  for (const auto i : make_range(_number_slip_systems)) {
	  
    _rho_ssd_increment[i] = _k_0 - _y_c * _rho_ssd[_qp][i];
    _rho_ssd_increment[i] *= std::abs(_slip_increment[_qp][i]) / _burgers_vector_mag;

  }
}

// calculate the irradiation dislocation loops increment based on equation (21)
void
CrystalPlasticityFerriticSteel::calculateDLincrement()
{

  // note that _slip_increment here is the rate
  for (const auto i : make_range(_number_slip_systems)) {

    _C_DL_increment[i] = (-1.0) * _lambda_DL * std::abs(_slip_increment[_qp][i])
                       * _C_DL[_qp][i] * _C_DL_diameter / _burgers_vector_mag;

  }
}

// calculate the irradiation solute cluster increment based on equation (23)
void
CrystalPlasticityFerriticSteel::calculateSCincrement() {

  // note that _slip_increment here is the rate
  for (const auto i : make_range(_number_slip_systems)) {

    _C_SC_increment[i] = (-1.0) * _lambda_SC * std::abs(_slip_increment[_qp][i])
                       * _C_SC[_qp][i];

  }
}

bool
CrystalPlasticityFerriticSteel::updateStateVariables()
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

  // irradiation dislocation loops
  for (const auto i : make_range(_number_slip_systems))
  {
    _C_DL_increment[i] *= _substep_dt;
    _C_DL[_qp][i] = _previous_substep_C_DL[i] + _C_DL_increment[i];
  }

  // solute clusters induced by irradiation
  for (const auto i : make_range(_number_slip_systems))
  {
    _C_SC_increment[i] *= _substep_dt;
    _C_SC[_qp][i] = _previous_substep_C_SC[i] + _C_SC_increment[i];
  }

  return true;
}
