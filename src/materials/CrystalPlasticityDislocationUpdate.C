// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// University of Bristol
// 27 Marzo 2022

#include "CrystalPlasticityDislocationUpdate.h"
#include "libmesh/int_range.h"
#include <cmath>

registerMooseObject("TensorMechanicsApp", CrystalPlasticityDislocationUpdate);

InputParameters
CrystalPlasticityDislocationUpdate::validParams()
{
  InputParameters params = CrystalPlasticityDislocationUpdateBase::validParams();
  params.addClassDescription("Dislocation based model for crystal plasticity "
                             "using the stress update code");
  params.addParam<Real>("ao", 0.001, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");  
  params.addParam<Real>("burgers_vector_mag",0.000256,"Magnitude of the Burgers vector");
  params.addParam<Real>("shear_modulus",86000.0,"Shear modulus in Taylor hardening law G");
  params.addParam<Real>("alpha_0",0.3,"Prefactor of Taylor hardening law, alpha");
  params.addParam<Real>("r", 1.4, "Latent hardening coefficient");
  params.addParam<Real>("tau_c_0", 0.112, "Peierls stress");
  params.addParam<Real>("k_0",100.0,"Coefficient K in SSD evolution, representing accumulation rate");
  params.addParam<Real>("y_c",0.0026,"Critical annihilation diameter");
  params.addParam<Real>("rho_tol",1.0,"Tolerance on dislocation density update");
  params.addParam<Real>("init_rho_ssd",1.0,"Initial dislocation density");
  params.addParam<Real>("init_rho_gnd_edge",0.0,"Initial dislocation density");
  params.addParam<Real>("init_rho_gnd_screw",0.0,"Initial dislocation density");
  params.addParam<MaterialPropertyName>(
      "total_twin_volume_fraction",
      "Total twin volume fraction, if twinning is considered in the simulation");
  return params;
}

CrystalPlasticityDislocationUpdate::CrystalPlasticityDislocationUpdate(
    const InputParameters & parameters)
  : CrystalPlasticityDislocationUpdateBase(parameters),
  
    // Constitutive model parameters
    _ao(getParam<Real>("ao")),
    _xm(getParam<Real>("xm")),
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")),
	_shear_modulus(getParam<Real>("shear_modulus")),
	_alpha_0(getParam<Real>("alpha_0")),
    _r(getParam<Real>("r")),
	_tau_c_0(getParam<Real>("tau_c_0")),
	_k_0(getParam<Real>("k_0")),
	_y_c(getParam<Real>("y_c")),
	
	// Initial values of the state variables
    _init_rho_ssd(getParam<Real>("init_rho_ssd")),
    _init_rho_gnd_edge(getParam<Real>("init_rho_gnd_edge")),
    _init_rho_gnd_screw(getParam<Real>("init_rho_gnd_screw")),
	
	// Tolerance on dislocation density update
	_rho_tol(getParam<Real>("rho_tol")),
	
	// not needed anymore, we only need slip resistance
    //_slip_resistance_increment(_number_slip_systems, 0.0),
	
    _rho_ssd(declareProperty<std::vector<Real>>("rho_ssd")),
    _rho_ssd_old(getMaterialPropertyOld<std::vector<Real>>("rho_ssd")),
    _rho_gnd_edge(declareProperty<std::vector<Real>>("rho_gnd_edge")),
   	_rho_gnd_edge_old(getMaterialPropertyOld<std::vector<Real>>("rho_gnd_edge")),
  	_rho_gnd_screw(declareProperty<std::vector<Real>>("rho_gnd_screw")),
    _rho_gnd_screw_old(getMaterialPropertyOld<std::vector<Real>>("rho_gnd_screw")),

    // resize local caching vectors used for substepping
    //_previous_substep_slip_resistance(_number_slip_systems, 0.0),
    //_slip_resistance_before_update(_number_slip_systems, 0.0),
	
    _rho_ssd_increment(_number_slip_systems, 0.0),
    _rho_gnd_edge_increment(_number_slip_systems, 0.0),
    _rho_gnd_screw_increment(_number_slip_systems, 0.0),
	
	// resize local caching vectors used for substepping
    _previous_substep_rho_ssd(_number_slip_systems, 0.0),
	_previous_substep_rho_gnd_edge(_number_slip_systems, 0.0),
	_previous_substep_rho_gnd_screw(_number_slip_systems, 0.0),
    _rho_ssd_before_update(_number_slip_systems, 0.0),
    _rho_gnd_edge_before_update(_number_slip_systems, 0.0),
    _rho_gnd_screw_before_update(_number_slip_systems, 0.0),  	

    // Twinning contributions, if used
    _include_twinning_in_Lp(parameters.isParamValid("total_twin_volume_fraction")),
    _twin_volume_fraction_total(_include_twinning_in_Lp
                                    ? &getMaterialPropertyOld<Real>("total_twin_volume_fraction")
                                    : nullptr),

    // store edge and screw slip directions to calculate directional derivatives
    // of the plastic slip rate	
    _edge_slip_direction(declareProperty<std::vector<Real>>("edge_slip_direction")),
	_screw_slip_direction(declareProperty<std::vector<Real>>("screw_slip_direction"))
{
}

void
CrystalPlasticityDislocationUpdate::initQpStatefulProperties()
{
  // Slip resistance is resized here
  CrystalPlasticityDislocationUpdateBase::initQpStatefulProperties();
  
  Real taylor_hardening;

  // Initialize and dislocation density size
  _rho_ssd[_qp].resize(_number_slip_systems);
  _rho_gnd_edge[_qp].resize(_number_slip_systems);
  _rho_gnd_screw[_qp].resize(_number_slip_systems);
  
  // Initialize dislocation densities
  for (const auto i : make_range(_number_slip_systems))
  {
    _rho_ssd[_qp][i] = _init_rho_ssd;
    _rho_gnd_edge[_qp][i] = _init_rho_gnd_edge;
    _rho_gnd_screw[_qp][i] = _init_rho_gnd_screw;
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

// Calculate Schmid tensor and
// store edge and screw slip directions to calculate directional derivatives
// of the plastic slip rate
void
CrystalPlasticityDislocationUpdate::calculateSchmidTensor(
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
CrystalPlasticityDislocationUpdate::setInitialConstitutiveVariableValues()
{
  // Would also set old dislocation densities here if included in this model
  //_slip_resistance[_qp] = _slip_resistance_old[_qp];
  _rho_ssd[_qp] = _rho_ssd_old[_qp];
  _previous_substep_rho_ssd = _rho_ssd_old[_qp];
  _rho_gnd_edge[_qp] = _rho_gnd_edge_old[_qp];
  _previous_substep_rho_gnd_edge = _rho_gnd_edge_old[_qp];
  _rho_gnd_screw[_qp] = _rho_gnd_screw_old[_qp];
  _previous_substep_rho_gnd_screw = _rho_gnd_screw_old[_qp];

  //_previous_substep_slip_resistance = _slip_resistance_old[_qp];
}

void
CrystalPlasticityDislocationUpdate::setSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  //_slip_resistance[_qp] = _previous_substep_slip_resistance;
  _rho_ssd[_qp] = _previous_substep_rho_ssd;
  _rho_gnd_edge[_qp] = _previous_substep_rho_gnd_edge;
  _rho_gnd_screw[_qp] = _previous_substep_rho_gnd_screw;
}

// Slip resistance can be calculated from dislocation density here only
// because it is the first method in which it is used,
// while calculateConstitutiveSlipDerivative is called afterwards
bool
CrystalPlasticityDislocationUpdate::calculateSlipRate()
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

// Slip resistance based on Taylor hardening
void
CrystalPlasticityDislocationUpdate::calculateSlipResistance()
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
	
	//std::cout << _slip_resistance[_qp][0] << std::endl;
	
	
}

void
CrystalPlasticityDislocationUpdate::calculateEquivalentSlipIncrement(
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
CrystalPlasticityDislocationUpdate::calculateConstitutiveSlipDerivative(
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
CrystalPlasticityDislocationUpdate::areConstitutiveStateVariablesConverged()
{
  return isConstitutiveStateVariableConverged(_rho_ssd[_qp],
                                              _rho_ssd_before_update,
                                              _previous_substep_rho_ssd,
                                              _rho_tol);

  // How do we check the tolerance of GNDs and is it needed?
											  
}

void
CrystalPlasticityDislocationUpdate::updateSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _previous_substep_rho_ssd = _rho_ssd[_qp];
  _previous_substep_rho_gnd_edge = _rho_gnd_edge[_qp];
  _previous_substep_rho_gnd_screw = _rho_gnd_screw[_qp];
}

void
CrystalPlasticityDislocationUpdate::cacheStateVariablesBeforeUpdate()
{
  _rho_ssd_before_update = _rho_ssd[_qp];
  _rho_gnd_edge_before_update = _rho_gnd_edge[_qp];
  _rho_gnd_screw_before_update = _rho_gnd_screw[_qp];
}

void
CrystalPlasticityDislocationUpdate::calculateStateVariableEvolutionRateComponent()
{
  Real rho_sum;

  // SSD dislocation density increment
  for (const auto i : make_range(_number_slip_systems))
  {
    
    rho_sum = _rho_ssd[_qp][i] + std::abs(_rho_gnd_edge[_qp][i]) + std::abs(_rho_gnd_screw[_qp][i]);

    // Multiplication and annihilation
	// note that _slip_increment here is the rate
	// and is multiplied by time step in updateStateVariables
    _rho_ssd_increment[i] = _k_0 * sqrt(rho_sum) - 2 * _y_c * _rho_ssd[_qp][i];
    _rho_ssd_increment[i] *= std::abs(_slip_increment[_qp][i]) / _burgers_vector_mag;

  }
  
  // GND dislocation density increment
  for (const auto i : make_range(_number_slip_systems)) {
	  
    // TO DO
    _rho_gnd_edge_increment[i] = 0.0;
    _rho_gnd_screw_increment[i] = 0.0;

    // no need for this, they are already assigned 
    // calculateSchmidTensor();
 
    // There is no need here for _edge_slip_direction
	// _grad_slip_increment is already the directional derivative
	// meaning \nabla \gamma \dot \hat{t} 
	// this is the rate
    //_rho_gnd_edge_increment[j] = (-1) * _grad_slip_increment[i] *  / _burgers_vector - _rho_gnd_edge_old[j];
    //_rho_gnd_screw_increment[j] = (-1) * _grad_slip_increment[i] * _screw_slip_direction[_qp][i * LIBMESH_DIM + j] / _burgers_vector - _rho_gnd_screw_old[j];    

  }

}

bool
CrystalPlasticityDislocationUpdate::updateStateVariables()
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
