// Nicolò Grilli
// University of Bristol
// Daijun Hu
// National University of Singapore
// 2 Agosto 2021

#include "DGAdvectionCoupledPN.h"

registerMooseObject("MooseApp", DGAdvectionCoupledPN);

InputParameters
DGAdvectionCoupledPN::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription("DG upwinding for the advection of a coupled variable. "
                             "Upwind condition is calculated both on edge/screw dislocations "
							 "in this element and on the neighbouring element. "
							 "The forward and backward motion of "
							 "positive and negative GND is taken into account.");
  params.addCoupledVar("rho_coupled", 0.0, "Coupled dislocation density in the flux term.");
  params.addCoupledVar("rho_coupled_ot", 0.0, "Other type: screw for edge kernel and vice versa.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  params.addParam<bool>("is_edge_or_screw",false,"Check if is edge or screw, or total density");
  return params;
}

DGAdvectionCoupledPN::DGAdvectionCoupledPN(const InputParameters & parameters)
  : DGKernel(parameters),
    _rho_coupled(coupledValue("rho_coupled")), // Coupled dislocation density in the flux term 
    _rho_coupled_coupled(isCoupled("rho_coupled")),
    _rho_coupled_var(_rho_coupled_coupled ? coupled("rho_coupled") : 0),
	_rho_coupled_neighbor(coupledNeighborValue("rho_coupled")),
    _rho_coupled_ot(coupledValue("rho_coupled_ot")), // Other type: screw for edge kernel and vice versa 
    _rho_coupled_ot_coupled(isCoupled("rho_coupled_ot")),
    _rho_coupled_ot_var(_rho_coupled_ot_coupled ? coupled("rho_coupled_ot") : 0),
	_rho_coupled_ot_neighbor(coupledNeighborValue("rho_coupled_ot")),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
    _dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")), // Velocity value (signed)
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>()),
    _is_edge_or_screw(getParam<bool>("is_edge_or_screw"))	
{
}

// read dislocation velocity from material object
// and store in _velocity
void
DGAdvectionCoupledPN::getDislocationVelocity()
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
  
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  {
	_velocity(j) *= _dislo_velocity[_qp][_slip_sys_index]; // velocity value (signed)
  }
	
}

Real
DGAdvectionCoupledPN::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0; // output residual
  Real rho_gnd_pos; // positive GND density in this element
  Real rho_gnd_neg; // negative GND density in this element
  Real neigh_rho_gnd_pos; // positive GND density in the neighbouring element
  Real neigh_rho_gnd_neg; // negative GND density in the neighbouring element
  
  // Remaining total dislocation density once the "other" gnd type is subtracted
  // This is necessary to consider the case in which both edge and screw GND are present
  Real remain_rho_tot;
  Real remain_rho_tot_neigh;
  
  getDislocationVelocity();

  Real vdotn = _velocity * _normals[_qp];
  
  // Define positive and negative GND densities
  // both are positive quantities
  
  if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
  
    // _rho_coupled[_qp] is rho_tot in this case
	// rho_t = sqrt(rho_e^2 + rho_s^2) in the pure GND case
	
	remain_rho_tot = _rho_coupled[_qp]*_rho_coupled[_qp]
	               - _rho_coupled_ot[_qp]*_rho_coupled_ot[_qp];

	if (remain_rho_tot >= 0.0) { 
		
      rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot) + _u[_qp]);
				  
	  rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot) - _u[_qp]);

	  rho_gnd_pos = std::max(rho_gnd_pos,0.0);
      rho_gnd_neg = std::max(rho_gnd_neg,0.0);

	} else { // All GNDs are of the other type, nothing left for this type

	  rho_gnd_pos = 0.0;
	  rho_gnd_neg = 0.0;

	}

	remain_rho_tot_neigh = _rho_coupled_neighbor[_qp]*_rho_coupled_neighbor[_qp]
	                     - _rho_coupled_ot_neighbor[_qp]*_rho_coupled_ot_neighbor[_qp];
						 
    if (remain_rho_tot_neigh >= 0.0) {
		
	  neigh_rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot_neigh) 
					    + _u_neighbor[_qp]);

	  neigh_rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot_neigh)  
	                    - _u_neighbor[_qp]);

      neigh_rho_gnd_pos = std::max(neigh_rho_gnd_pos,0.0);
	  neigh_rho_gnd_neg = std::max(neigh_rho_gnd_neg,0.0);
		
	} else {
		
	  neigh_rho_gnd_pos = 0.0;
	  neigh_rho_gnd_neg = 0.0;

	}
 
  } else // Case with rho_edge or rho_screw = _rho_coupled[_qp]
  {
	// _u[_qp] is rho_tot in this case
	// rho_t = sqrt(rho_e^2 + rho_s^2) in the pure GND case
	
	remain_rho_tot = _u[_qp]*_u[_qp] 
	               - _rho_coupled_ot[_qp]*_rho_coupled_ot[_qp];
	
	if (remain_rho_tot >= 0.0) {

      rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot) + _rho_coupled[_qp]);
      rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot) - _rho_coupled[_qp]);

	  rho_gnd_pos = std::max(rho_gnd_pos,0.0);
      rho_gnd_neg = std::max(rho_gnd_neg,0.0);	  
		
	} else {

	  rho_gnd_pos = 0.0;
	  rho_gnd_neg = 0.0;
		
	}
	
    remain_rho_tot_neigh = _u_neighbor[_qp]*_u_neighbor[_qp]
	                     - _rho_coupled_ot_neighbor[_qp]*_rho_coupled_ot_neighbor[_qp];
	
	if (remain_rho_tot_neigh >= 0.0) {
		
      neigh_rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot_neigh) 
						+ _rho_coupled_neighbor[_qp]);
						
      neigh_rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot_neigh) 
	                    - _rho_coupled_neighbor[_qp]);

      neigh_rho_gnd_pos = std::max(neigh_rho_gnd_pos,0.0);
	  neigh_rho_gnd_neg = std::max(neigh_rho_gnd_neg,0.0);  

	} else {
		
      neigh_rho_gnd_pos = 0.0;
	  neigh_rho_gnd_neg = 0.0;
		
	}
	  
  }

  switch (type)
  {
    case Moose::Element:	
	  
	  if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	    if (vdotn >= 0.0) {
          r += vdotn * rho_gnd_pos * _test[_i][_qp]; // positive gnd exits from element
		  r += vdotn * neigh_rho_gnd_neg * _test[_i][_qp]; // negative gnd enters from neighbour
		}
	  
	    if (vdotn < 0.0) {
		  r += vdotn * neigh_rho_gnd_pos * _test[_i][_qp]; // positive gnd enters from neighbour
		  r += vdotn * rho_gnd_neg * _test[_i][_qp]; // negative gnd exits from element
		}
		  
	  } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
		  
	    if (vdotn >= 0.0) {
          r += vdotn * rho_gnd_pos * _test[_i][_qp]; // rho_total exits from element
		  r -= vdotn * neigh_rho_gnd_neg * _test[_i][_qp]; // rho_total enters from neighbour
		}		  
	  
	    if (vdotn < 0.0) {
		  r += vdotn * neigh_rho_gnd_pos * _test[_i][_qp]; // rho_total enters from neighbour
	      r -= vdotn * rho_gnd_neg * _test[_i][_qp]; // rho_total exits from element
		}		  
	  }

      break;

    case Moose::Neighbor: // opposite signs than Moose::Element and use _test_neighbor[_i][_qp]
	
	  if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	    if (vdotn >= 0.0) {
          r -= vdotn * rho_gnd_pos * _test_neighbor[_i][_qp]; // positive gnd exits from element
		  r -= vdotn * neigh_rho_gnd_neg * _test_neighbor[_i][_qp]; // negative gnd enters from neighbour
		}
	  
	    if (vdotn < 0.0) {
		  r -= vdotn * neigh_rho_gnd_pos * _test_neighbor[_i][_qp]; // positive gnd enters from neighbour
		  r -= vdotn * rho_gnd_neg * _test_neighbor[_i][_qp]; // negative gnd exits from element
		}
		  
	  } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
	  
	  
	    if (vdotn >= 0.0) {
          r -= vdotn * rho_gnd_pos * _test_neighbor[_i][_qp]; // rho_total exits from element
		  r += vdotn * neigh_rho_gnd_neg * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
		}		  
	  
	    if (vdotn < 0.0) {
		  r -= vdotn * neigh_rho_gnd_pos * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
	      r += vdotn * rho_gnd_neg * _test_neighbor[_i][_qp]; // rho_total exits from element
		}
	  }
	  
      break;
  }

  return r;
}

// Jacobian is non-zero because both _u and rho_coupled are in the residual
Real 
DGAdvectionCoupledPN::computeQpJacobian(Moose::DGJacobianType type) 
{ 
  Real r = 0; // output Jacobian
  
  getDislocationVelocity();

  Real vdotn = _velocity * _normals[_qp];
  
  switch (type)
  {
    case Moose::ElementElement:	
	  
	  if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	    if (vdotn >= 0.0) {
          r += vdotn * 0.5 * _phi[_j][_qp] * _test[_i][_qp]; // positive gnd exits from element
		}
	  
	    if (vdotn < 0.0) {
		  r += vdotn * (-0.5) * _phi[_j][_qp] * _test[_i][_qp]; // negative gnd exits from element
		}
		  
	  } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
		  
	    if (vdotn >= 0.0) {
          r += vdotn * 0.5 * _phi[_j][_qp] * _test[_i][_qp]; // rho_total exits from element
		}		  
	  
	    if (vdotn < 0.0) {
	      r -= vdotn * 0.5 * _phi[_j][_qp] * _test[_i][_qp]; // rho_total exits from element
		}		  
	  }

      break;	  
	  
    case Moose::ElementNeighbor:	
	  
	  if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	    if (vdotn >= 0.0) {
		  r += vdotn * (-0.5) * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // negative gnd enters from neighbour
		}
	  
	    if (vdotn < 0.0) {
		  r += vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // positive gnd enters from neighbour
		}
		  
	  } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
		  
	    if (vdotn >= 0.0) {
		  r -= vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // rho_total enters from neighbour
		}		  
	  
	    if (vdotn < 0.0) {
		  r += vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // rho_total enters from neighbour
		}		  
	  }

      break;

    case Moose::NeighborElement:
	
	  if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	    if (vdotn >= 0.0) {
          r -= vdotn * 0.5 * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // positive gnd exits from element
		}
	  
	    if (vdotn < 0.0) {
		  r -= vdotn * (-0.5) * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // negative gnd exits from element
		}
		  
	  } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
	  
	  
	    if (vdotn >= 0.0) {
          r -= vdotn * 0.5 * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total exits from element
		}		  
	  
	    if (vdotn < 0.0) {
	      r += vdotn * 0.5 * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total exits from element
		}
	  }
	  
      break;
	    
	case Moose::NeighborNeighbor:
	
	  if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	    if (vdotn >= 0.0) {
		  r -= vdotn * (-0.5) * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // negative gnd enters from neighbour
		}
	  
	    if (vdotn < 0.0) {
		  r -= vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // positive gnd enters from neighbour
		}
		  
	  } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
	  
	  
	    if (vdotn >= 0.0) {
		  r += vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
		}		  
	  
	    if (vdotn < 0.0) {
		  r -= vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
		}
	  }	
	
	  break;  
  }  

  return r; 
}

Real
DGAdvectionCoupledPN::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real r = 0; // output off diagonal Jacobian
  Real vdotn;

  if (_rho_coupled_coupled && jvar == _rho_coupled_var) {
	  
    getDislocationVelocity();

    vdotn = _velocity * _normals[_qp];

    switch (type)
    {
      case Moose::ElementElement:
	  
	    if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	      if (vdotn >= 0.0) {
            r += vdotn * 0.5 * _phi[_j][_qp] * _test[_i][_qp]; // positive gnd exits from element
		  }
	  
	      if (vdotn < 0.0) {
		    r += vdotn * 0.5 * _phi[_j][_qp] * _test[_i][_qp]; // negative gnd exits from element
		  }
		  
	    } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
		  
	      if (vdotn >= 0.0) {
            r += vdotn * 0.5 * _phi[_j][_qp] * _test[_i][_qp]; // rho_total exits from element
		  }		  
	  
	      if (vdotn < 0.0) {
	        r -= vdotn * (-0.5) * _phi[_j][_qp] * _test[_i][_qp]; // rho_total exits from element
		  }		  
	    }
		
        break;

      case Moose::ElementNeighbor:
	  
	    if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	      if (vdotn >= 0.0) {
		    r += vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // negative gnd enters from neighbour
		  }
	  
	      if (vdotn < 0.0) {
		    r += vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // positive gnd enters from neighbour
		  }
		  
	    } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
		  
	      if (vdotn >= 0.0) {
		    r -= vdotn * (-0.5) * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // rho_total enters from neighbour
		  }		  
	  
	      if (vdotn < 0.0) {
		    r += vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test[_i][_qp]; // rho_total enters from neighbour
		  }		  
	    }
		
        break;

      case Moose::NeighborElement:
	  
	    if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	      if (vdotn >= 0.0) {
            r -= vdotn * 0.5 * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // positive gnd exits from element
		  }
	  
	      if (vdotn < 0.0) {
		    r -= vdotn * 0.5 * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // negative gnd exits from element
		  }
		  
	    } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
	  
	      if (vdotn >= 0.0) {
            r -= vdotn * 0.5 * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total exits from element
		  }		  
	  
	      if (vdotn < 0.0) {
	        r += vdotn * (-0.5) * _phi[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total exits from element
		  }
	    }
		
        break;

      case Moose::NeighborNeighbor:
	  
	    if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	      if (vdotn >= 0.0) {
		    r -= vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // negative gnd enters from neighbour
		  }
	  
	      if (vdotn < 0.0) {
		    r -= vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // positive gnd enters from neighbour
		  }
		  
	    } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
	  
	      if (vdotn >= 0.0) {
		    r += vdotn * (-0.5) * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
		  }		  
	  
	      if (vdotn < 0.0) {
		    r -= vdotn * 0.5 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp]; // rho_total enters from neighbour
		  }
	    }
		
        break;

    }

  }

  return r;
}
