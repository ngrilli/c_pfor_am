// Nicolo Grilli
// University of Bristol
// Daijun Hu
// National University of Singapore
// 20 Ottobre 2021

#include "DGA2Trace.h"

registerMooseObject("MooseApp", DGA2Trace);

InputParameters
DGA2Trace::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription("Dislocation curvature diffusion "
                             "second term on the right hand side of equation (3) in "
							 "Stefan Sandfeld and Michael Zaiser "
							 "Pattern formation in a minimal model of continuum dislocation plasticity "
							 "Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp). "
							 "Discontinuous Galerkin formulation");
  params.addCoupledVar("rho_tot", 0.0, "Total dislocation density: rho_t.");
  params.addCoupledVar("dv_dx", 0.0, "Derivative of the velocity with respect to edge slip direction.");
  params.addCoupledVar("dv_dy", 0.0, "Derivative of the velocity with respect to screw slip direction.");
  params.addParam<Real>("dv_dx_max",1e9,"Max absolute value of dv_dx");
  params.addParam<Real>("dv_dy_max",1e9,"Max absolute value of dv_dy");    
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");  
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  return params;
}

DGA2Trace::DGA2Trace(const InputParameters & parameters)
  : DGKernel(parameters),
    _rho_tot(coupledValue("rho_tot")), // Total dislocation density: rho_t
    _rho_tot_coupled(isCoupled("rho_tot")),
    _rho_tot_var(_rho_tot_coupled ? coupled("rho_tot") : 0),
	_rho_neighbor(coupledNeighborValue("rho_tot")),
    _dv_dx(coupledValue("dv_dx")), // Derivative of the velocity with respect to edge slip direction
    _dv_dy(coupledValue("dv_dy")), // Derivative of the velocity with respect to screw slip direction
	_dv_dx_neighbor(coupledNeighborValue("dv_dx")), // same in the neighbouring element
	_dv_dy_neighbor(coupledNeighborValue("dv_dy")), // same in the neighbouring element
    _dv_dx_max(getParam<Real>("dv_dx_max")),
	_dv_dy_max(getParam<Real>("dv_dy_max")),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>())
{
}

// read dislocation velocity direction from material object
// and store in _velocity
// note that multiplication by velocity magnitude is not carried out
// because the term inside the derivatives d/dx and d/dy
// does not include the velocity magnitude
void
DGA2Trace::getDislocationVelocity()
{
	
  // Find dislocation velocity based on slip systems index and dislocation character
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	  {
	    _velocity(j) = _edge_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // edge direction	  
	  }
	  break;
	case DisloCharacter::screw:
	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	  {
		// note that the definition of _screw_slip_direction in FiniteStrainCrystalPlasticityDislo
		// is -y, because +x is _edge_slip_direction and +z is slip plane normal
		// but derivative must be taken along +y
		// therefore a sign change is needed
	    _velocity(j) = - _screw_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // screw direction	  
	  }	
	  break;
  }
	
}

Real
DGA2Trace::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0.0;
  Real advected_quantity; // term inside the d/dx or d/dy derivative
  Real neigh_advected_quantity; // same in the neighbouring element
  
  // Velocity derivatives
  Real dv_dx;
  Real dv_dy;
  Real neigh_dv_dx;
  Real neigh_dv_dy;
  
  getDislocationVelocity();

  Real vdotn = _velocity * _normals[_qp];
  
  // Limit value of velocity derivatives
  dv_dx = std::min(std::abs(_dv_dx[_qp]),_dv_dx_max);
  dv_dy = std::min(std::abs(_dv_dy[_qp]),_dv_dy_max);
  
  dv_dx = dv_dx * std::copysign(1.0, _dv_dx[_qp]);
  dv_dy = dv_dy * std::copysign(1.0, _dv_dy[_qp]);
  
  // Same in the neighbouring element
  neigh_dv_dx = std::min(std::abs(_dv_dx_neighbor[_qp]),_dv_dx_max);
  neigh_dv_dy = std::min(std::abs(_dv_dy_neighbor[_qp]),_dv_dy_max);
  
  neigh_dv_dx = neigh_dv_dx * std::copysign(1.0, _dv_dx_neighbor[_qp]);
  neigh_dv_dy = neigh_dv_dy * std::copysign(1.0, _dv_dy_neighbor[_qp]);
  
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
      advected_quantity = 0.5 * _rho_tot[_qp] * dv_dx;
	  neigh_advected_quantity = 0.5 * _rho_neighbor[_qp] * neigh_dv_dx;
	  break;
	case DisloCharacter::screw:
      advected_quantity = 0.5 * _rho_tot[_qp] * dv_dy;
	  neigh_advected_quantity = 0.5 * _rho_neighbor[_qp] * neigh_dv_dy;
	  break;	  
  }

  switch (type)
  {
    case Moose::Element:
      if (vdotn >= 0)
        r += vdotn * advected_quantity * _test[_i][_qp]; // flux out of this element
      else
        r += vdotn * neigh_advected_quantity * _test[_i][_qp]; // flux in from neighbour
      break;

    case Moose::Neighbor: // opposite sign in neighbour residual to obtain conservative flux
      if (vdotn >= 0)
        r -= vdotn * advected_quantity * _test_neighbor[_i][_qp]; // flux out of this element
      else
        r -= vdotn * neigh_advected_quantity * _test_neighbor[_i][_qp]; // flux in from neighbour
      break;
  }

  return r;
}

// This term does not depend on the total curvature q_t
// therefore in-diagonal Jacobian is zero
Real
DGA2Trace::computeQpJacobian(Moose::DGJacobianType type)
{
  return 0.0;
}

// This term depend on the total dislocation density rho_t
// So calculate derivative with respect to rho_t
Real
DGA2Trace::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real r = 0.0;
  Real d_advected; // d advected_quantity /d rho_tot
  Real d_neigh_advected; // // d neigh_advected_quantity /d rho_neighbor
  
  // Velocity derivatives
  Real dv_dx;
  Real dv_dy;
  Real neigh_dv_dx;
  Real neigh_dv_dy;
  
  // scalar product between direction (edge or screw) 
  // and interface normal between this element and the neighbour
  Real vdotn;
  
  if (_rho_tot_coupled && jvar == _rho_tot_var)
  {
  
    getDislocationVelocity();

    vdotn = _velocity * _normals[_qp];
  
    // Limit value of velocity derivatives
    dv_dx = std::min(std::abs(_dv_dx[_qp]),_dv_dx_max);
    dv_dy = std::min(std::abs(_dv_dy[_qp]),_dv_dy_max);
  
    dv_dx = dv_dx * std::copysign(1.0, _dv_dx[_qp]);
    dv_dy = dv_dy * std::copysign(1.0, _dv_dy[_qp]);
  
    // Same in the neighbouring element
    neigh_dv_dx = std::min(std::abs(_dv_dx_neighbor[_qp]),_dv_dx_max);
    neigh_dv_dy = std::min(std::abs(_dv_dy_neighbor[_qp]),_dv_dy_max);
  
    neigh_dv_dx = neigh_dv_dx * std::copysign(1.0, _dv_dx_neighbor[_qp]);
    neigh_dv_dy = neigh_dv_dy * std::copysign(1.0, _dv_dy_neighbor[_qp]);
  
    switch (_dislo_character)
    {
      case DisloCharacter::edge:
        d_advected = 0.5 * _phi[_j][_qp] * dv_dx;
	    d_neigh_advected = 0.5 * _phi_neighbor[_j][_qp] * neigh_dv_dx;
	    break;
	  case DisloCharacter::screw:
        d_advected = 0.5 * _phi[_j][_qp] * dv_dy;
	    d_neigh_advected = 0.5 * _phi_neighbor[_j][_qp] * neigh_dv_dy;
	    break;	  
    }

    switch (type)
    {
      case Moose::ElementElement:
	  
        if (vdotn >= 0.0) {
          r += vdotn * d_advected * _test[_i][_qp];
        }
        break;

      case Moose::ElementNeighbor:

        if (vdotn < 0) {
          r += vdotn * d_neigh_advected * _test[_i][_qp];
        }
        break;

      case Moose::NeighborElement:
      
        if (vdotn >= 0) {
          r -= vdotn * d_advected * _test_neighbor[_i][_qp];
        }
        break;
	    
      case Moose::NeighborNeighbor:
         
		if (vdotn < 0) {
          r -= vdotn * d_neigh_advected * _test_neighbor[_i][_qp];
        }
        break;
	 
    } // end of switch case
  
  } // end of if

  return r;
  
}
