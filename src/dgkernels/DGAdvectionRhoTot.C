// Nicolò Grilli
// University of Bristol
// Daijun Hu
// National University of Singapore
// 4 Settembre 2021

#include "DGAdvectionRhoTot.h"

registerMooseObject("MooseApp", DGAdvectionRhoTot);

InputParameters
DGAdvectionRhoTot::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription("DG upwinding for the advection of a coupled variable. "
                             "Upwind condition is calculated both on edge/screw dislocations "
							 "in this element and on the neighbouring element. "
							 "The forward and backward motion of "
							 "positive and negative GND is taken into account. "
							 "The velocity directions of GND with arbitrary character "
							 "perpendicular to the line segment is taken into account. "
							 "This kernel must be applied to rho_tot.");
  params.addCoupledVar("rho_edge", 0.0, "Edge dislocation density in the flux term.");
  params.addCoupledVar("rho_screw", 0.0, "Screw dislocation density in the flux term.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw. "
									 "This kernel implements: "
									 "d(rho_gnd v)/dx if dislo_character = edge "
									 "d(rho_gnd v)/dy if dislo_character = screw.");
  params.addParam<bool>("check_gnd_rho_ratio",false,"Check that |rho_gnd| / rho_tot <= 1");
  params.addParam<Real>("rho_tot_tol",0.000001,"Tolerance on small values of rho_tot.");
  return params;
}

DGAdvectionRhoTot::DGAdvectionRhoTot(const InputParameters & parameters)
  : DGKernel(parameters),
    _rho_edge(coupledValue("rho_edge")), // Coupled dislocation density in the flux term 
    _rho_edge_coupled(isCoupled("rho_edge")),
    _rho_edge_var(_rho_edge_coupled ? coupled("rho_edge") : 0),
	_rho_edge_neighbor(coupledNeighborValue("rho_edge")),
    _rho_screw(coupledValue("rho_screw")), // Other type: screw for edge kernel and vice versa 
    _rho_screw_coupled(isCoupled("rho_screw")),
    _rho_screw_var(_rho_screw_coupled ? coupled("rho_screw") : 0),
	_rho_screw_neighbor(coupledNeighborValue("rho_screw")),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
    _dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")), // Velocity value (signed)
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>()),
	_check_gnd_rho_ratio(getParam<bool>("check_gnd_rho_ratio")),
    _rho_tot_tol(getParam<Real>("rho_tot_tol"))	// Tolerance on small values of rho_tot
{
}

// read dislocation velocity from material object
// and store in _velocity
void
DGAdvectionRhoTot::getDislocationVelocity()
{
  RealVectorValue edge_velocity;
  RealVectorValue screw_velocity;
  
  // Ratios between GND densities and total dislocation density in this element
  Real edge_rho_tot_ratio = 0.0;
  Real screw_rho_tot_ratio = 0.0;
  
  // Ratios between GND densities and total dislocation density in the neighbouring element
  Real edge_rho_tot_ratio_neigh = 0.0;
  Real screw_rho_tot_ratio_neigh = 0.0;
  
  // Final angle cosines used to calculate velocity direction
  // theta = 0 means pure edge GND
  Real costheta;
  Real sintheta;
  Real onetheta;
  
  // Allocate dislocation velocities based on slip systems index and dislocation character
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j) 
  {
    edge_velocity(j) = _edge_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // edge direction	  	
    screw_velocity(j) = - _screw_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // screw direction
  }
  
  // Find ratio between GND densities and total dislocation density in this element
  if (_u[_qp] > _rho_tot_tol) {
	  
	edge_rho_tot_ratio = _rho_edge[_qp] / _u[_qp];
	screw_rho_tot_ratio = _rho_screw[_qp] / _u[_qp];

    if (_check_gnd_rho_ratio) {
		
	  if (std::abs(edge_rho_tot_ratio) > 1.0) {
		
        edge_rho_tot_ratio = std::copysign(1.0, _rho_edge[_qp]);		
			
      }

      if (std::abs(screw_rho_tot_ratio) > 1.0) {
		  
	    screw_rho_tot_ratio = std::copysign(1.0, _rho_screw[_qp]);
		  
	  }	  
		
	}	  
	  
  }

  // Find ratio between GND densities and total dislocation density in the neighbouring element
  if (_u_neighbor[_qp] > _rho_tot_tol) {
	  
	edge_rho_tot_ratio_neigh = _rho_edge_neighbor[_qp] / _u_neighbor[_qp];
	screw_rho_tot_ratio_neigh = _rho_screw_neighbor[_qp] / _u_neighbor[_qp];

    if (_check_gnd_rho_ratio) {
		
	  if (std::abs(edge_rho_tot_ratio_neigh) > 1.0) {
		
        edge_rho_tot_ratio_neigh = std::copysign(1.0, _rho_edge_neighbor[_qp]);		
			
      }

      if (std::abs(screw_rho_tot_ratio_neigh) > 1.0) {
		  
	    screw_rho_tot_ratio_neigh = std::copysign(1.0, _rho_screw_neighbor[_qp]);
		  
	  }	  
		
	}	  
	  
  }
  
  // Ensure that the magnitude of the dislocation velocity
  // depends only on the load in case
  // edge_rho_tot_ratio and screw_rho_tot_ratio are too low
  // If GND densities in the current element are too low
  // check the values from the neighbouring element
  if (std::abs(edge_rho_tot_ratio) > _rho_tot_tol || 
      std::abs(screw_rho_tot_ratio) > _rho_tot_tol) 
  {
	  
	onetheta = std::sqrt(edge_rho_tot_ratio * edge_rho_tot_ratio + 
	                     screw_rho_tot_ratio * screw_rho_tot_ratio);
	  
    costheta = edge_rho_tot_ratio / onetheta;
	sintheta = screw_rho_tot_ratio / onetheta;				
	  
  } else if (std::abs(edge_rho_tot_ratio_neigh) > _rho_tot_tol || 
             std::abs(screw_rho_tot_ratio_neigh) > _rho_tot_tol)
  {
	
	onetheta = std::sqrt(edge_rho_tot_ratio_neigh * edge_rho_tot_ratio_neigh + 
	                     screw_rho_tot_ratio_neigh * screw_rho_tot_ratio_neigh);
	  
    costheta = edge_rho_tot_ratio_neigh / onetheta;
	sintheta = screw_rho_tot_ratio_neigh / onetheta;
	
  } else 
  {
	  
	// Otherwise assume it is a GND dislocation
	// with character dislo_character
    switch (_dislo_character)
    {
      case DisloCharacter::edge:
	
        costheta = 1.0;
	    sintheta = 0.0;

	    break;
	  
	  case DisloCharacter::screw:
	
        costheta = 0.0;
	    sintheta = 1.0;

	    break;
    }	
	
  }

  // Find dislocation velocity based on GND state
  // The multiplication by velocity value (signed)
  // ensures that dislocation direction is consistent with the load
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  {
    _velocity(j) = costheta * edge_velocity(j);
	_velocity(j) += sintheta * screw_velocity(j);
  }  
  
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  {
	_velocity(j) *= _dislo_velocity[_qp][_slip_sys_index]; // velocity value (signed)
  }
  
  // Velocity needs to be multiplied by costheta or sintheta
  // depending on character because in case of pure GND
  // the two kernels DGAdvectionRhoTot for edge and screw character
  // must give a total contribution to the flux proportional to 
  // rho_t costheta^2 + rho_t sintheta^2
  // and not to
  // rho_t costheta + rho_t sintheta 
  // vdotn does not account for that factor because it is a sum of
  // edge and screw velocity directions weighted with costheta and sintheta  
  // This is not needed in DGAdvectionRhoGND because the quantity advected is rho_gnd itself
  // and there is only one kernel
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
	
	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      {
        _velocity(j) *= costheta;
	  }

	  break;
	  
	case DisloCharacter::screw:

	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      {	
        _velocity(j) *= sintheta;
	  }

	  break;
  }

}

Real
DGAdvectionRhoTot::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0; // output residual
  
  // Edge or screw depending on dislo_character
  Real rho_gnd_pos; // positive GND density in this element
  Real rho_gnd_neg; // negative GND density in this element
  Real neigh_rho_gnd_pos; // positive edge GND density in the neighbouring element
  Real neigh_rho_gnd_neg; // negative edge GND density in the neighbouring element
  
  // Remaining total density after subtraction of "other type" GND
  Real remain_rho_tot;
  Real remain_rho_tot_neigh;
  
  // Coupled GND density:
  // _rho_edge[_qp] if dislo_character = edge
  // _rho_screw[_qp] if dislo_character = screw
  Real rho_coupled;
  Real rho_coupled_neigh;
  
  getDislocationVelocity();

  Real vdotn = _velocity * _normals[_qp];
  
  // Define positive and negative GND densities
  // both are positive quantities
  // edge or screw depending on dislo_character
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
	
      // rho_t = sqrt(rho_e^2 + rho_s^2) in the pure GND case
      remain_rho_tot = _u[_qp]*_u[_qp] - _rho_screw[_qp]*_rho_screw[_qp];
	  
	  rho_coupled = _rho_edge[_qp];
	  
	  remain_rho_tot_neigh = _u_neighbor[_qp]*_u_neighbor[_qp]
	                       - _rho_screw_neighbor[_qp]*_rho_screw_neighbor[_qp];
						   
	  rho_coupled_neigh = _rho_edge_neighbor[_qp];

	  break;
	  
	case DisloCharacter::screw:
	
      // rho_t = sqrt(rho_e^2 + rho_s^2) in the pure GND case
      remain_rho_tot = _u[_qp]*_u[_qp] - _rho_edge[_qp]*_rho_edge[_qp];
	  
	  rho_coupled = _rho_screw[_qp];
	  
	  remain_rho_tot_neigh = _u_neighbor[_qp]*_u_neighbor[_qp]
	                       - _rho_edge_neighbor[_qp]*_rho_edge_neighbor[_qp];
						   
	  rho_coupled_neigh = _rho_screw_neighbor[_qp];

	  break;
  }

  if (remain_rho_tot >= 0.0) { 
		
    rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot) + rho_coupled);
				 
    rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot) - rho_coupled);

    rho_gnd_pos = std::max(rho_gnd_pos,0.0);
    rho_gnd_neg = std::max(rho_gnd_neg,0.0);

  } else { // All GNDs are of the other type, nothing left for this type

    rho_gnd_pos = 0.0;
	rho_gnd_neg = 0.0;

  }
						 
  if (remain_rho_tot_neigh >= 0.0) {
		
    neigh_rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot_neigh) 
					  + rho_coupled_neigh);

    neigh_rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot_neigh)  
	                  - rho_coupled_neigh);

    neigh_rho_gnd_pos = std::max(neigh_rho_gnd_pos,0.0);
	neigh_rho_gnd_neg = std::max(neigh_rho_gnd_neg,0.0);
		
  } else {
		
	neigh_rho_gnd_pos = 0.0;
	neigh_rho_gnd_neg = 0.0;

  }

  // Calculate residual
  switch (type)
  {
    case Moose::Element:
	
	  if (vdotn >= 0.0) {
        r += vdotn * rho_gnd_pos * _test[_i][_qp]; // rho_total exits from element
		r -= vdotn * neigh_rho_gnd_neg * _test[_i][_qp]; // rho_total enters from neighbour
      }		  
	  
	  if (vdotn < 0.0) {
        r += vdotn * neigh_rho_gnd_pos * _test[_i][_qp]; // rho_total enters from neighbour
	    r -= vdotn * rho_gnd_neg * _test[_i][_qp]; // rho_total exits from element
      }		  

      break;

    case Moose::Neighbor: // opposite signs than Moose::Element and use _test_neighbor[_i][_qp]

	  if (vdotn >= 0.0) {
        r -= vdotn * rho_gnd_pos * _test_neighbor[_i][_qp]; // rho_total exits from element
		r += vdotn * neigh_rho_gnd_neg * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
	  }		  
	  
	  if (vdotn < 0.0) {
		r -= vdotn * neigh_rho_gnd_pos * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
	    r += vdotn * rho_gnd_neg * _test_neighbor[_i][_qp]; // rho_total exits from element
      }
	  
      break;
  }

  return r;
}

// Jacobian is non-zero because both _u and rho_coupled are in the residual
Real 
DGAdvectionRhoTot::computeQpJacobian(Moose::DGJacobianType type) 
{ 
  Real r = 0; // output Jacobian
  
  getDislocationVelocity();

  Real vdotn = _velocity * _normals[_qp];
  
  switch (type)
  {
    case Moose::ElementElement:
		  
	  if (vdotn >= 0.0) {
        r += vdotn * 0.5 * _phi[_j][_qp] * _test[_i][_qp]; // rho_total exits from element
	  }		  
	  
	  if (vdotn < 0.0) {
	    r -= vdotn * 0.5 * _phi[_j][_qp] * _test[_i][_qp]; // rho_total exits from element
	  }

      break;	  
	  
    case Moose::ElementNeighbor:	
		  
	  if (vdotn >= 0.0) {
	    r -= vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // rho_total enters from neighbour
	  }		  
	  
	  if (vdotn < 0.0) {
        r += vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // rho_total enters from neighbour
	  }

      break;

    case Moose::NeighborElement:
	  
	  if (vdotn >= 0.0) {
        r -= vdotn * 0.5 * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total exits from element
	  }		  
	  
	  if (vdotn < 0.0) {
	    r += vdotn * 0.5 * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total exits from element
	  }
	  
      break;
	    
	case Moose::NeighborNeighbor:

	  if (vdotn >= 0.0) {
	    r += vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
	  }		  
	  
	  if (vdotn < 0.0) {
        r -= vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
	  }	
	
	  break;  
  }  

  return r; 
}

// Derivatives with respect to _rho_edge[_qp] and _rho_screw[_qp] 
Real
DGAdvectionRhoTot::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real r = 0; // output off diagonal Jacobian
  Real vdotn;

  // The following equations are valid  
  // for both d/dx kernel and d/dy kernel
  if ((_rho_edge_coupled && jvar == _rho_edge_var) || 
      (_rho_screw_coupled && jvar == _rho_screw_var)) {
		  
    getDislocationVelocity();

    vdotn = _velocity * _normals[_qp];

    switch (type)
    {
      case Moose::ElementElement:
	  
	    if (vdotn >= 0.0) {
          r += vdotn * 0.5 * _phi[_j][_qp] * _test[_i][_qp]; // rho_total exits from element
		}		  
	  
	    if (vdotn < 0.0) {
	      r -= vdotn * (-0.5) * _phi[_j][_qp] * _test[_i][_qp]; // rho_total exits from element
		}

        break;

      case Moose::ElementNeighbor:
		  
	    if (vdotn >= 0.0) {
		  r -= vdotn * (-0.5) * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // rho_total enters from neighbour
		}		  
	  
	    if (vdotn < 0.0) {
		  r += vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // rho_total enters from neighbour
		}
		
        break;

      case Moose::NeighborElement:
	  
	    if (vdotn >= 0.0) {
          r -= vdotn * 0.5 * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total exits from element
		}		  
	  
	    if (vdotn < 0.0) {
	      r += vdotn * (-0.5) * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total exits from element
		}
		
        break;

      case Moose::NeighborNeighbor:
	  
	    if (vdotn >= 0.0) {
		  r += vdotn * (-0.5) * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
		}		  
	  
	    if (vdotn < 0.0) {
		  r -= vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
		}
		
        break;

    }

  }

  return r;
}
