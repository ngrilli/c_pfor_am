// Nicolo Grilli
// National University of Singapore
// 25 Gennaio 2021

#include "CurvatureAdvection.h"

registerMooseObject("MooseApp", CurvatureAdvection);

defineLegacyParams(CurvatureAdvection);

InputParameters
CurvatureAdvection::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Dislocation curvature advection "
                             "first term on the right hand side of equation (3) in "
							 "Stefan Sandfeld and Michael Zaiser "
							 "Pattern formation in a minimal model of continuum dislocation plasticity "
							 "Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp). ");
  params.addCoupledVar("rho_gnd", 0.0, "GND dislocation density: rho_x or rho_y for edge or screw.");
  params.addCoupledVar("rho_tot", 0.0, "Total dislocation density: rho_t.");
  params.addParam<Real>("rho_tot_tol",0.000001,"Tolerance on small values of rho_tot.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  MooseEnum dislo_sign("positive negative", "positive");
  params.addRequiredParam<MooseEnum>("dislo_sign",
                                     dislo_sign,
                                     "Sign of dislocations.");
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  return params;
}

CurvatureAdvection::CurvatureAdvection(const InputParameters & parameters)
  : Kernel(parameters),
    _rho_gnd(coupledValue("rho_gnd")), // GND dislocation density: rho_x or rho_y for edge or screw
    _rho_gnd_coupled(isCoupled("rho_gnd")),
    _rho_gnd_var(_rho_gnd_coupled ? coupled("rho_gnd") : 0),
    _rho_tot(coupledValue("rho_tot")), // Total dislocation density: rho_t
    _rho_tot_coupled(isCoupled("rho_tot")),
    _rho_tot_var(_rho_tot_coupled ? coupled("rho_tot") : 0),
	_rho_tot_tol(getParam<Real>("rho_tot_tol")),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
    _dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")), // Velocity value (signed)
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dislo_sign(getParam<MooseEnum>("dislo_sign").getEnum<DisloSign>()),
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>())
{
}

Real
CurvatureAdvection::negSpeedQp()
{
  Real edge_sign;

  switch (_dislo_sign)
  {
    case DisloSign::positive:
      edge_sign = 1.0;
      break;
    case DisloSign::negative:
      edge_sign = -1.0;
      break;
  }  
	
  _velocity.resize(3, 0.0);
  
  // Find dislocation velocity based on slip systems index and dislocation character
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	  {
	    _velocity[j] = _edge_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // edge direction	  
	  }
	  break;
	case DisloCharacter::screw:
	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	  {
		// note that the definition of _screw_slip_direction in FiniteStrainCrystalPlasticityDislo
		// is -y, because +x is _edge_slip_direction and +z is slip plane normal
		// but derivative must be taken along +y
		// therefore a sign change is needed
	    _velocity[j] = - _screw_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // screw direction	  
	  }	
	  break;
  }
	
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  {
	_velocity[j] *= _dislo_velocity[_qp][_slip_sys_index]; // velocity value (signed)
	_velocity[j] *= edge_sign; // positive or negative dislocation
  }
  
  return -_grad_test[_i][_qp] * RealVectorValue(_velocity[0],_velocity[1],_velocity[2]);
}

// This is the no-upwinded version
Real
CurvatureAdvection::computeQpResidual()
{
  Real val;
  
  if (_rho_tot[_qp] > _rho_tot_tol) {
	  
	val = _rho_gnd[_qp] * _u[_qp] / _rho_tot[_qp];  
	
  } else {
	  
	val = 0.0;  
	  
  }

  return negSpeedQp() * val;
}

// This is the no-upwinded version
Real
CurvatureAdvection::computeQpJacobian()
{
  Real val;
  
  if (_rho_tot[_qp] > _rho_tot_tol) {
	  
    val = _rho_gnd[_qp] / _rho_tot[_qp];  
	  
  } else {
	  
    val = 0.0;  	  
	  
  }
  
  return negSpeedQp() * _phi[_j][_qp] * val;
}

// This is the no-upwinded version
Real
CurvatureAdvection::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_rho_gnd_coupled && jvar == _rho_gnd_var)
  {
	
	if (_rho_tot[_qp] > _rho_tot_tol) {
	
      return negSpeedQp() * _phi[_j][_qp] * _u[_qp] / _rho_tot[_qp];
	
	} else {
		
	  return 0.0;	
		
	}
	
  } else if (_rho_tot_coupled && jvar == _rho_tot_var) {
	  
	if (_rho_tot[_qp] > _rho_tot_tol) {
	  
	  return (-1.0) * negSpeedQp() * _phi[_j][_qp] * _rho_gnd[_qp] * _u[_qp] / (_rho_tot[_qp] * _rho_tot[_qp]);
	
	} else {
		
	  return 0.0;	
		
	}
	
  } else {
	  
	return 0.0;  
	  
  }
}
