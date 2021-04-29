// Nicolo Grilli
// National University of Singapore
// 27 Gennaio 2021

#include "A2Trace.h"

registerMooseObject("MooseApp", A2Trace);

defineLegacyParams(A2Trace);

InputParameters
A2Trace::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Dislocation curvature diffusion "
                             "second term on the right hand side of equation (3) in "
							 "Stefan Sandfeld and Michael Zaiser "
							 "Pattern formation in a minimal model of continuum dislocation plasticity "
							 "Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp). ");
  params.addCoupledVar("rho_tot", 0.0, "Total dislocation density: rho_t.");
  params.addCoupledVar("dv_dx", 0.0, "Derivative of the velocity with respect to edge slip direction.");
  params.addCoupledVar("dv_dy", 0.0, "Derivative of the velocity with respect to screw slip direction.");
  params.addParam<Real>("dv_dx_max",1e9,"Max absolute value of dv_dx");
  params.addParam<Real>("dv_dy_max",1e9,"Max absolute value of dv_dy");  
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

A2Trace::A2Trace(const InputParameters & parameters)
  : Kernel(parameters),
    _rho_tot(coupledValue("rho_tot")), // Total dislocation density: rho_t
    _rho_tot_coupled(isCoupled("rho_tot")),
    _rho_tot_var(_rho_tot_coupled ? coupled("rho_tot") : 0),
    _dv_dx(coupledValue("dv_dx")), // Derivative of the velocity with respect to edge slip direction
    _dv_dy(coupledValue("dv_dy")), // Derivative of the velocity with respect to screw slip direction
    _dv_dx_max(getParam<Real>("dv_dx_max")),
	_dv_dy_max(getParam<Real>("dv_dy_max")),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dislo_sign(getParam<MooseEnum>("dislo_sign").getEnum<DisloSign>()),
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>())
{
}

Real
A2Trace::negSpeedQp()
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
	_velocity[j] *= edge_sign; // positive or negative dislocation
  }
  
  return -_grad_test[_i][_qp] * RealVectorValue(_velocity[0],_velocity[1],_velocity[2]);
}

// This is the no-upwinded version
Real
A2Trace::computeQpResidual()
{
  Real val;
  Real dv_dx;
  Real dv_dy;
  
  dv_dx = std::min(std::abs(_dv_dx[_qp]),_dv_dx_max);
  dv_dy = std::min(std::abs(_dv_dy[_qp]),_dv_dy_max);
  
  dv_dx = dv_dx * std::copysign(1.0, _dv_dx[_qp]);
  dv_dy = dv_dy * std::copysign(1.0, _dv_dy[_qp]);
  
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
      val = 0.5 * _rho_tot[_qp] * dv_dx;
	  break;
	case DisloCharacter::screw:
      val = 0.5 * _rho_tot[_qp] * dv_dy;
	  break;	  
  }

  return negSpeedQp() * val;
}

// This is the no-upwinded version
Real
A2Trace::computeQpJacobian()
{
  // no dependence on the curvature
  return 0.0;
}

// This is the no-upwinded version
Real
A2Trace::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  Real dv_dx;
  Real dv_dy;
  
  dv_dx = std::min(std::abs(_dv_dx[_qp]),_dv_dx_max);
  dv_dy = std::min(std::abs(_dv_dy[_qp]),_dv_dy_max);
  
  dv_dx = dv_dx * std::copysign(1.0, _dv_dx[_qp]);
  dv_dy = dv_dy * std::copysign(1.0, _dv_dy[_qp]);
  
  if (_rho_tot_coupled && jvar == _rho_tot_var)
  {
	
    switch (_dislo_character)
    {
      case DisloCharacter::edge:
        val = 0.5 * _phi[_j][_qp] * dv_dx;
	    break;
	  case DisloCharacter::screw:
        val = 0.5 * _phi[_j][_qp] * dv_dy;
	    break;	  
    }
	
    return negSpeedQp() * val;

  } else {
	  
	return 0.0;  
	  
  }
}
