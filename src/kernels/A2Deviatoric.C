// Nicolo Grilli
// National University of Singapore
// 30 Gennaio 2021

#include "A2Deviatoric.h"

#include <cmath>

registerMooseObject("MooseApp", A2Deviatoric);

defineLegacyParams(A2Deviatoric);

InputParameters
A2Deviatoric::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Dislocation curvature diffusion "
                             "second term on the right hand side of equation (3) in "
							 "Stefan Sandfeld and Michael Zaiser (deviatoric part) "
							 "Pattern formation in a minimal model of continuum dislocation plasticity "
							 "Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp). ");
  params.addCoupledVar("rho_gnd_edge", 0.0, "Edge dislocation density: rho_x.");
  params.addCoupledVar("rho_gnd_screw", 0.0, "Screw dislocation density: rho_y.");  
  params.addCoupledVar("dv_dx", 0.0, "Derivative of the velocity with respect to edge slip direction.");
  params.addCoupledVar("dv_dy", 0.0, "Derivative of the velocity with respect to screw slip direction.");
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

A2Deviatoric::A2Deviatoric(const InputParameters & parameters)
  : Kernel(parameters),
    _rho_gnd_edge(coupledValue("rho_gnd_edge")), // Edge dislocation density: rho_x
    _rho_gnd_edge_coupled(isCoupled("rho_gnd_edge")),
    _rho_gnd_edge_var(_rho_gnd_edge_coupled ? coupled("rho_gnd_edge") : 0),
    _rho_gnd_screw(coupledValue("rho_gnd_screw")), // Screw dislocation density: rho_y
    _rho_gnd_screw_coupled(isCoupled("rho_gnd_screw")),
    _rho_gnd_screw_var(_rho_gnd_screw_coupled ? coupled("rho_gnd_screw") : 0),	
    _dv_dx(coupledValue("dv_dx")), // Derivative of the velocity with respect to edge slip direction
    _dv_dy(coupledValue("dv_dy")), // Derivative of the velocity with respect to screw slip direction	
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dislo_sign(getParam<MooseEnum>("dislo_sign").getEnum<DisloSign>()),
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>())
{
}

Real
A2Deviatoric::negSpeedQp()
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
A2Deviatoric::computeQpResidual()
{
  Real val; // return value
  Real ksabs; // modulus of k GND vector
  Real diagterm; // diagonal term of deviatoric A2 matrix
  Real outofdiagterm; // otu of diagonal term of deviatoric A2 matrix
  
  // sqrt(rho_x^2 + rho_y^2)
  ksabs = std::sqrt(_rho_gnd_edge[_qp]*_rho_gnd_edge[_qp]+_rho_gnd_screw[_qp]*_rho_gnd_screw[_qp]);
  
  // rho_x^2 - rho_y^2
  diagterm = _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp] - _rho_gnd_screw[_qp]*_rho_gnd_screw[_qp];
  
  // rho_x rho_y
  outofdiagterm = _rho_gnd_edge[_qp]*_rho_gnd_screw[_qp];
  
  if (ksabs > 0.001) {
    diagterm = 0.5 * diagterm / ksabs;
	outofdiagterm = outofdiagterm / ksabs;
  } else {
	diagterm = 0.0;  
	outofdiagterm = 0.0;
  }
  
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
	  val = (-1.0) * diagterm * _dv_dx[_qp];
	  val = val - outofdiagterm * _dv_dy[_qp];
	  break;
	case DisloCharacter::screw:
      val = diagterm * _dv_dy[_qp];
	  val = val - outofdiagterm * _dv_dx[_qp];
	  break;	  
  }

  return negSpeedQp() * val;
}

// This is the no-upwinded version
Real
A2Deviatoric::computeQpJacobian()
{
  // no dependence on the curvature
  return 0.0;
}

// This is the no-upwinded version
Real
A2Deviatoric::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  Real ksabs3; // modulus of k GND vector, power 3  
  Real ddiagterm_drhox;
  Real doutofdiagterm_drhox;
  Real ddiagterm_drhoy;
  Real doutofdiagterm_drhoy;
  
  // (rho_x^2 + rho_y^2)^(3/2)
  ksabs3 = std::pow(_rho_gnd_edge[_qp]*_rho_gnd_edge[_qp]+_rho_gnd_screw[_qp]*_rho_gnd_screw[_qp],1.5);
  
  if (_rho_gnd_edge_coupled && jvar == _rho_gnd_edge_var)
  {
	// derivative with respect to rho_x
	
	ddiagterm_drhox = _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp]+3.0*_rho_gnd_screw[_qp]*_rho_gnd_screw[_qp];
	ddiagterm_drhox = 0.5 * _rho_gnd_edge[_qp] * ddiagterm_drhox;
	
	if (ksabs3 > 0.001) {
	  ddiagterm_drhox = ddiagterm_drhox / ksabs3;
	} else {
	  ddiagterm_drhox = 0.0;	
	}
		
	doutofdiagterm_drhox = _rho_gnd_screw[_qp]*_rho_gnd_screw[_qp]*_rho_gnd_screw[_qp];
	
	if (ksabs3 > 0.001) {
	  doutofdiagterm_drhox = doutofdiagterm_drhox / ksabs3;
	} else {
	  doutofdiagterm_drhox = 0.0;	
	}
	
    switch (_dislo_character)
    {
      case DisloCharacter::edge:
	  
	    // Residual is:
	  	// - diagterm * _dv_dx[_qp] - outofdiagterm * _dv_dy[_qp]

        val = (-1.0) * ddiagterm_drhox * _dv_dx[_qp] - doutofdiagterm_drhox * _dv_dy[_qp];
		
	    break;
		
	  case DisloCharacter::screw:
	  
	  	// Residual is:
	  	// - outofdiagterm * _dv_dx[_qp] + diagterm * _dv_dy[_qp]
		
        val = (-1.0) * doutofdiagterm_drhox * _dv_dx[_qp] + ddiagterm_drhox * _dv_dy[_qp];
		
	    break;	  
    }
	
    return negSpeedQp() * _phi[_j][_qp] * val;

  } else if (_rho_gnd_screw_coupled && jvar == _rho_gnd_screw_var) {
	  
	// derivative with respect to rho_y
	
	ddiagterm_drhoy = _rho_gnd_screw[_qp]*_rho_gnd_screw[_qp]+3.0*_rho_gnd_edge[_qp]*_rho_gnd_edge[_qp];
	ddiagterm_drhoy = (-0.5) * _rho_gnd_screw[_qp] * ddiagterm_drhoy;
	
	if (ksabs3 > 0.001) {
	  ddiagterm_drhoy = ddiagterm_drhoy / ksabs3;
	} else {
      ddiagterm_drhoy = 0.0;
	}
		
	doutofdiagterm_drhoy = _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp]*_rho_gnd_edge[_qp];
	
	if (ksabs3 > 0.001) {
	  doutofdiagterm_drhoy = doutofdiagterm_drhoy / ksabs3;
	} else {
      doutofdiagterm_drhoy = 0.0;
	}	
	
    switch (_dislo_character)
    {
      case DisloCharacter::edge:
	  
	    // Residual is:
	  	// - diagterm * _dv_dx[_qp] - outofdiagterm * _dv_dy[_qp]
		
        val = (-1.0) * ddiagterm_drhoy * _dv_dx[_qp] - doutofdiagterm_drhoy * _dv_dy[_qp];
		
	    break;
		
	  case DisloCharacter::screw:
	  
	  	// Residual is:
	  	// - outofdiagterm * _dv_dx[_qp] + diagterm * _dv_dy[_qp]
	  
        val = (-1.0) * doutofdiagterm_drhoy * _dv_dx[_qp] + ddiagterm_drhoy * _dv_dy[_qp];
		
	    break;	  
    }
	
    return negSpeedQp() * _phi[_j][_qp] * val;	  
	  
  } else {
	  
	return 0.0;  
	  
  }
}
