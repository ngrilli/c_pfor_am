// Nicolo Grilli
// University of Bristol
// Daijun Hu
// National University of Singapore
// 4 Agosto 2021

#include "DGCurvatureAdvectionPN.h"

registerMooseObject("MooseApp", DGCurvatureAdvectionPN);

InputParameters
DGCurvatureAdvectionPN::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription("DG Dislocation curvature advection."
							 "Upwind condition is calculated both on edge/screw dislocations "
							 "in this element and on the neighbouring element. "
							 "The forward and backward motion of "
							 "positive and negative GND is taken into account. ");
  params.addCoupledVar("rho_gnd", 0.0, "GND dislocation density: rho_x or rho_y for edge or screw.");
  params.addCoupledVar("rho_gnd_ot", 0.0, "Other type: screw for edge kernel and vice versa.");
  params.addCoupledVar("rho_tot", 0.0, "Total dislocation density: rho_t.");
  params.addParam<Real>("rho_tot_tol",0.000001,"Tolerance on small values of rho_tot.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  params.addParam<bool>("check_gnd_rho_ratio",false,"Check that |rho_gnd| / rho_tot <= 1");
  return params;
}

DGCurvatureAdvectionPN::DGCurvatureAdvectionPN(const InputParameters & parameters)
  : DGKernel(parameters),
    _rho_gnd(coupledValue("rho_gnd")), // GND dislocation density: rho_x or rho_y for edge or screw
    _rho_gnd_coupled(isCoupled("rho_gnd")),
    _rho_gnd_var(_rho_gnd_coupled ? coupled("rho_gnd") : 0),
    _rho_tot(coupledValue("rho_tot")), // Total dislocation density: rho_t
    _rho_tot_coupled(isCoupled("rho_tot")),
    _rho_tot_var(_rho_tot_coupled ? coupled("rho_tot") : 0),
	_rho_gnd_neighbor(coupledNeighborValue("rho_gnd")),
    _rho_tot_neighbor(coupledNeighborValue("rho_tot")),
    _rho_gnd_ot(coupledValue("rho_gnd_ot")), // Other type: screw for edge kernel and vice versa 
    _rho_gnd_ot_coupled(isCoupled("rho_gnd_ot")),
    _rho_gnd_ot_var(_rho_gnd_ot_coupled ? coupled("rho_gnd_ot") : 0),
	_rho_gnd_ot_neighbor(coupledNeighborValue("rho_gnd_ot")),
	_rho_tot_tol(getParam<Real>("rho_tot_tol")), // tolerance on small values of the total dislocation density
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
    _dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")), // Velocity value (signed)
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>()),
	_check_gnd_rho_ratio(getParam<bool>("check_gnd_rho_ratio"))
{
}

// read dislocation velocity from material object
// and store in _velocity
void
DGCurvatureAdvectionPN::getDislocationVelocity()
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
DGCurvatureAdvectionPN::computeQpResidual(Moose::DGResidualType type)
{ 
  Real r = 0; // Residual for output
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
  // _u[_qp] is q_tot in this case
  // _rho_tot[_qp] is rho_tot in this case
  // _rho_gnd[_qp] is edge or screw, respectively
  // _rho_gnd_ot[_qp] is screw or edge, respectively
  // rho_t = sqrt(rho_e^2 + rho_s^2) in the pure GND case
  
  remain_rho_tot = _rho_tot[_qp]*_rho_tot[_qp]
	               - _rho_gnd_ot[_qp]*_rho_gnd_ot[_qp];
				   
  if (remain_rho_tot >= 0.0) {
	  
    rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot) + _rho_gnd[_qp]);
				  
    rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot) - _rho_gnd[_qp]);

    rho_gnd_pos = std::max(rho_gnd_pos,0.0);
    rho_gnd_neg = std::max(rho_gnd_neg,0.0);		  
 
  } else { // All GNDs are of the other type, nothing left for this type
	  
	rho_gnd_pos = 0.0;
	rho_gnd_neg = 0.0;
	  
  }
  
  remain_rho_tot_neigh = _rho_tot_neighbor[_qp]*_rho_tot_neighbor[_qp]
	                     - _rho_gnd_ot_neighbor[_qp]*_rho_gnd_ot_neighbor[_qp];
				 
  if (remain_rho_tot_neigh >= 0.0) {
	  
    neigh_rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot_neigh) + _rho_gnd_neighbor[_qp]);
						
    neigh_rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot_neigh) - _rho_gnd_neighbor[_qp]);

    neigh_rho_gnd_pos = std::max(neigh_rho_gnd_pos,0.0);
    neigh_rho_gnd_neg = std::max(neigh_rho_gnd_neg,0.0);

  } else { // All GNDs are of the other type, nothing left for this type
		
    neigh_rho_gnd_pos = 0.0;
    neigh_rho_gnd_neg = 0.0;
	
  }						 
  
  // residual calculation: this is the curvature residual
  // if a positive curvature exits the element, residual must be positive 
  // if a positive curvature enters the element, residual must be negative
  // There is no dependence on the sign of the gnd considered  

  switch (type)
  {
    case Moose::Element:
	  	  
	  if (vdotn >= 0.0) {
		  
		if (_rho_tot[_qp] > _rho_tot_tol) {
			
          // positive gnd exits from element
          r += ((vdotn * rho_gnd_pos * _u[_qp]) / _rho_tot[_qp]) 
		       * _test[_i][_qp];
			   
		}
		 
        if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
		
          // negative gnd enters from neighbour 		
		  r -= ((vdotn * neigh_rho_gnd_neg * _u_neighbor[_qp]) / _rho_tot_neighbor[_qp]) 
		       * _test[_i][_qp];
		
		}
		
	  }

      if (vdotn < 0.0) {
		  
        if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
			
		  // positive gnd enters from neighbour
	      r += ((vdotn * neigh_rho_gnd_pos * _u_neighbor[_qp]) / _rho_tot_neighbor[_qp]) 
	           * _test[_i][_qp];		
			
		}

        if (_rho_tot[_qp] > _rho_tot_tol) {
			
          // negative gnd exits from element
          r -= ((vdotn * rho_gnd_neg * _u[_qp]) / _rho_tot[_qp]) 
	           * _test[_i][_qp];				
			
		}			 

	  }
		  
      break;

    case Moose::Neighbor:
	  		  
	  if (vdotn >= 0.0) {
		  
		if (_rho_tot[_qp] > _rho_tot_tol) {
			
          // positive gnd exits from element
          r -= ((vdotn * rho_gnd_pos * _u[_qp]) / _rho_tot[_qp]) 
		       * _test_neighbor[_i][_qp];
			   
		}
		 
        if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
		
          // negative gnd enters from neighbour		
		  r += ((vdotn * neigh_rho_gnd_neg * _u_neighbor[_qp]) / _rho_tot_neighbor[_qp]) 
		       * _test_neighbor[_i][_qp];
		
		}
		
	  }

      if (vdotn < 0.0) {
		  
        if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
			
		  // positive gnd enters from neighbour
	      r -= ((vdotn * neigh_rho_gnd_pos * _u_neighbor[_qp]) / _rho_tot_neighbor[_qp]) 
	           * _test_neighbor[_i][_qp];		
			
		}

        if (_rho_tot[_qp] > _rho_tot_tol) {
			
          // negative gnd exits from element
          r += ((vdotn * rho_gnd_neg * _u[_qp]) / _rho_tot[_qp]) 
	           * _test_neighbor[_i][_qp];			
			
		}			 

	  }
		  
      break;
  }

  return r;
}

// Exact form of the in-diagonal Jacobian considering the signed gnd densities
Real 
DGCurvatureAdvectionPN::computeQpJacobian(Moose::DGJacobianType type) 
{ 
  Real jac = 0; // Jacobian for output
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
  // _u[_qp] is q_tot in this case
  // _rho_tot[_qp] is rho_tot in this case
  // _rho_gnd[_qp] is edge or screw, respectively
  // _rho_gnd_ot[_qp] is screw or edge, respectively
  // rho_t = sqrt(rho_e^2 + rho_s^2) in the pure GND case
  
  remain_rho_tot = _rho_tot[_qp]*_rho_tot[_qp]
	             - _rho_gnd_ot[_qp]*_rho_gnd_ot[_qp];
  
  if (remain_rho_tot >= 0.0) { 
  
    rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot) + _rho_gnd[_qp]);
				  
    rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot) - _rho_gnd[_qp]);

    rho_gnd_pos = std::max(rho_gnd_pos,0.0);
    rho_gnd_neg = std::max(rho_gnd_neg,0.0);	

  } else { // All GNDs are of the other type, nothing left for this type

	rho_gnd_pos = 0.0;
	rho_gnd_neg = 0.0;

  }	
  
  remain_rho_tot_neigh = _rho_tot_neighbor[_qp]*_rho_tot_neighbor[_qp]
	                   - _rho_gnd_ot_neighbor[_qp]*_rho_gnd_ot_neighbor[_qp];  

  if (remain_rho_tot_neigh >= 0.0) {
	  
    neigh_rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot_neigh) + _rho_gnd_neighbor[_qp]);
						
    neigh_rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot_neigh) - _rho_gnd_neighbor[_qp]);

    neigh_rho_gnd_pos = std::max(neigh_rho_gnd_pos,0.0);
    neigh_rho_gnd_neg = std::max(neigh_rho_gnd_neg,0.0);
	
  } else { // All GNDs are of the other type, nothing left for this type
		
	neigh_rho_gnd_pos = 0.0;
	neigh_rho_gnd_neg = 0.0;
	  
  }  
					
  switch (type)
  {
    case Moose::ElementElement:
	
      if (vdotn >= 0.0) {
		  
		if (_rho_tot[_qp] > _rho_tot_tol) {
			
          // positive gnd exits from element
          jac += ((vdotn * rho_gnd_pos * _phi[_j][_qp]) / _rho_tot[_qp]) 
		         * _test[_i][_qp];
			   
		}		  
		  
	  }	
	
	  if (vdotn < 0.0) {
		  
        if (_rho_tot[_qp] > _rho_tot_tol) {
			
          // negative gnd exits from element
          jac -= ((vdotn * rho_gnd_neg * _phi[_j][_qp]) / _rho_tot[_qp]) 
	             * _test[_i][_qp];				
			
		}			  
		  
	  }	
		  
      break;

    case Moose::ElementNeighbor:
	
	  if (vdotn >= 0.0) {
		  
        if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
		
          // negative gnd enters from neighbour 		
		  jac -= ((vdotn * neigh_rho_gnd_neg * _phi_neighbor[_j][_qp]) / _rho_tot_neighbor[_qp]) 
		         * _test[_i][_qp];
		
		}		  
		  
	  }
	
	  if (vdotn < 0.0) {
		  
        if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
			
		  // positive gnd enters from neighbour
	      jac += ((vdotn * neigh_rho_gnd_pos * _phi_neighbor[_j][_qp]) / _rho_tot_neighbor[_qp]) 
	             * _test[_i][_qp];		
			
		}		  
		  
	  }
		  
      break;
	  
	case Moose::NeighborElement:
	
	  if (vdotn >= 0.0) {
		  
		if (_rho_tot[_qp] > _rho_tot_tol) {
			
          // positive gnd exits from element
          jac -= ((vdotn * rho_gnd_pos * _phi[_j][_qp]) / _rho_tot[_qp]) 
		         * _test_neighbor[_i][_qp];
			   
		}		  
		  
	  }
	
	  if (vdotn < 0.0) {
		  
        if (_rho_tot[_qp] > _rho_tot_tol) {
			
          // negative gnd exits from element
          jac += ((vdotn * rho_gnd_neg * _phi[_j][_qp]) / _rho_tot[_qp]) 
	             * _test_neighbor[_i][_qp];			
			
		}		  
		  
	  }

	  break;
	
	case Moose::NeighborNeighbor:
	
      if (vdotn >= 0.0) {
		  
        if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
		
          // negative gnd enters from neighbour		
		  jac += ((vdotn * neigh_rho_gnd_neg * _phi_neighbor[_j][_qp]) / _rho_tot_neighbor[_qp]) 
		         * _test_neighbor[_i][_qp];
		
		}		  
		  
	  }	
	
	  if (vdotn < 0.0) {
		  
        if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
			
		  // positive gnd enters from neighbour
	      jac -= ((vdotn * neigh_rho_gnd_pos * _phi_neighbor[_j][_qp]) / _rho_tot_neighbor[_qp]) 
	             * _test_neighbor[_i][_qp];		
			
		}
		  
	  }
	
	  break;
  }  

  return jac; 
}

// An approximated form of the non-diagonal Jacobian
// with the approximation _rho_coupled_ot[_qp] = _rho_coupled_ot_neighbor[_qp] = 0
// so no derivative with respect to _rho_coupled_ot or _rho_coupled_ot_neighbor
Real
DGCurvatureAdvectionPN::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real jac = 0; // Jacobian for output
  Real rho_gnd_pos = 0; // positive GND density in this element
  Real rho_gnd_neg = 0; // negative GND density in this element
  Real neigh_rho_gnd_pos = 0; // positive GND density in the neighbouring element
  Real neigh_rho_gnd_neg = 0; // negative GND density in the neighbouring element
  
  // Remaining total dislocation density once the "other" gnd type is subtracted
  // This is necessary to consider the case in which both edge and screw GND are present  
  Real remain_rho_tot;
  Real remain_rho_tot_neigh;
  
  getDislocationVelocity();
  
  Real vdotn = _velocity * _normals[_qp];
  
  // Define positive and negative GND densities
  // both are positive quantities
  // _u[_qp] is q_tot in this case
  // _rho_tot[_qp] is rho_tot in this case
  // _rho_gnd[_qp] is edge or screw, respectively
  // _rho_gnd_ot[_qp] is screw or edge, respectively
  // rho_t = sqrt(rho_e^2 + rho_s^2) in the pure GND case
  
  if (_rho_tot_coupled && jvar == _rho_tot_var) {
	  
	// Positive and negative GND densities
    // are needed only for the derivative with
    // respect to _rho_tot[_qp]
  
    remain_rho_tot = _rho_tot[_qp]*_rho_tot[_qp]
	               - _rho_gnd_ot[_qp]*_rho_gnd_ot[_qp];
  
    if (remain_rho_tot >= 0.0) {
		
      rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot) + _rho_gnd[_qp]);
				  
      rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot) - _rho_gnd[_qp]);

      rho_gnd_pos = std::max(rho_gnd_pos,0.0);
      rho_gnd_neg = std::max(rho_gnd_neg,0.0);	

    } else { // All GNDs are of the other type, nothing left for this type

	  rho_gnd_pos = 0.0;
	  rho_gnd_neg = 0.0;

	}	
	
    remain_rho_tot_neigh = _rho_tot_neighbor[_qp]*_rho_tot_neighbor[_qp]
	                     - _rho_gnd_ot_neighbor[_qp]*_rho_gnd_ot_neighbor[_qp];  

    if (remain_rho_tot_neigh >= 0.0) {

      neigh_rho_gnd_pos = 0.5 * (std::sqrt(remain_rho_tot_neigh) + _rho_gnd_neighbor[_qp]);
						
      neigh_rho_gnd_neg = 0.5 * (std::sqrt(remain_rho_tot_neigh) - _rho_gnd_neighbor[_qp]);

      neigh_rho_gnd_pos = std::max(neigh_rho_gnd_pos,0.0);
      neigh_rho_gnd_neg = std::max(neigh_rho_gnd_neg,0.0);

    } else {
		
	  neigh_rho_gnd_pos = 0.0;
	  neigh_rho_gnd_neg = 0.0;
	  
    }
 
  }

  if (_rho_gnd_coupled && jvar == _rho_gnd_var) {
	  
    // derivative with respect to _rho_gnd and _rho_gnd_neighbor
	switch (type)
    {
      case Moose::ElementElement:
	  
        if (vdotn >= 0.0) {	
		
		  if (_rho_tot[_qp] > _rho_tot_tol) {
			
            // positive gnd exits from element
            jac += ((vdotn * 0.5 * _phi[_j][_qp] * _u[_qp]) / _rho_tot[_qp]) 
		         * _test[_i][_qp];
			   
		  }		
		
		}  
	  
	    if (vdotn < 0.0) {
		
          if (_rho_tot[_qp] > _rho_tot_tol) {
			
            // negative gnd exits from element
            jac -= ((vdotn * (-0.5) * _phi[_j][_qp] * _u[_qp]) / _rho_tot[_qp]) 
	             * _test[_i][_qp];				
			
		  }		
			
		}
		
        break;

      case Moose::ElementNeighbor:
	  
        if (vdotn >= 0.0) {
			 
          if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
		
            // negative gnd enters from neighbour 		
		    jac -= ((vdotn * (-0.5) * _phi_neighbor[_j][_qp] * _u_neighbor[_qp]) 
			       / _rho_tot_neighbor[_qp]) * _test[_i][_qp];
		
		  }	
			
		}		  
	  
	    if (vdotn < 0.0) {
			
          if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
			
		    // positive gnd enters from neighbour
	        jac += ((vdotn * 0.5 * _phi_neighbor[_j][_qp] * _u_neighbor[_qp]) 
			       / _rho_tot_neighbor[_qp]) * _test[_i][_qp];		
			
		  }			

		}  
		
        break;	  
	  
	  case Moose::NeighborElement:
	  
        if (vdotn >= 0.0) {
			
		  if (_rho_tot[_qp] > _rho_tot_tol) {
			
            // positive gnd exits from element
            jac -= ((vdotn * 0.5 * _phi[_j][_qp] * _u[_qp]) / _rho_tot[_qp]) 
		           * _test_neighbor[_i][_qp];
			   
		  }			

		}	  
	  
        if (vdotn < 0.0) {
			
          if (_rho_tot[_qp] > _rho_tot_tol) {
			
            // negative gnd exits from element
            jac += ((vdotn * (-0.5) * _phi[_j][_qp] * _u[_qp]) / _rho_tot[_qp]) 
	               * _test_neighbor[_i][_qp];			
			
		  }		
			
		}	  
		
        break;

	  case Moose::NeighborNeighbor:
	  
        if (vdotn >= 0.0) {
			
          if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
		
            // negative gnd enters from neighbour		
		    jac += ((vdotn * (-0.5) * _phi_neighbor[_j][_qp] * _u_neighbor[_qp]) 
			       / _rho_tot_neighbor[_qp]) * _test_neighbor[_i][_qp];
		
		  }			
			
		}	  
	  
	    if (vdotn < 0.0) {
			
          if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
			
		    // positive gnd enters from neighbour
	        jac -= ((vdotn * 0.5 * _phi_neighbor[_j][_qp] * _u_neighbor[_qp]) 
			       / _rho_tot_neighbor[_qp]) * _test_neighbor[_i][_qp];		
			
		  }			
			
		}	  
		
        break;	  
	}
  
  } else if (_rho_tot_coupled && jvar == _rho_tot_var) {
	  
    // derivative with respect to _rho_tot and _rho_tot_neighbor
	// the numerator has also a derivative with respect to 
	// _rho_tot and _rho_tot_neighbor but here it is neglected
	switch (type)
    {
	  case Moose::ElementElement:
	  
	    if (vdotn >= 0.0) {
		  
		  if (_rho_tot[_qp] > _rho_tot_tol) {
			
            // positive gnd exits from element
            jac += (-1.0) * ((vdotn * rho_gnd_pos * _u[_qp]) / (_rho_tot[_qp] * _rho_tot[_qp])) 
		           * _phi[_j][_qp] * _test[_i][_qp];
			   
		  }
		
	    }

        if (vdotn < 0.0) {

          if (_rho_tot[_qp] > _rho_tot_tol) {
			
            // negative gnd exits from element
            jac -= (-1.0) * ((vdotn * rho_gnd_neg * _u[_qp]) / (_rho_tot[_qp] * _rho_tot[_qp])) 
	               * _phi[_j][_qp] * _test[_i][_qp];				
			
		  }			 

	    }  
	  
	    break;
		
	  case Moose::ElementNeighbor:
	  	  
	    if (vdotn >= 0.0) {
		 
          if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
		
            // negative gnd enters from neighbour 		
		    jac -= (-1.0) * ((vdotn * neigh_rho_gnd_neg * _u_neighbor[_qp]) 
			       / (_rho_tot_neighbor[_qp] * _rho_tot_neighbor[_qp])) 
		           * _phi_neighbor[_j][_qp] * _test[_i][_qp];
		
		  }
		
	    }

        if (vdotn < 0.0) {
		  
          if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
			
		    // positive gnd enters from neighbour
	        jac += (-1.0) * ((vdotn * neigh_rho_gnd_pos * _u_neighbor[_qp]) 
			       / (_rho_tot_neighbor[_qp] * _rho_tot_neighbor[_qp])) 
	               * _phi_neighbor[_j][_qp] * _test[_i][_qp];		
			
		  }		 

	    }		  
	  
	    break;
		
	  case Moose::NeighborElement:
	  
	    if (vdotn >= 0.0) {
		  
		  if (_rho_tot[_qp] > _rho_tot_tol) {
			
            // positive gnd exits from element
            jac -= (-1.0) * ((vdotn * rho_gnd_pos * _u[_qp]) / (_rho_tot[_qp] * _rho_tot[_qp])) 
		           * _phi[_j][_qp] * _test_neighbor[_i][_qp];
			   
		  }
		
	    }

        if (vdotn < 0.0) {

          if (_rho_tot[_qp] > _rho_tot_tol) {
			
            // negative gnd exits from element
            jac += (-1.0) * ((vdotn * rho_gnd_neg * _u[_qp]) / (_rho_tot[_qp] * _rho_tot[_qp])) 
	             * _phi[_j][_qp] * _test_neighbor[_i][_qp];			
			
		  }			 

	    } 
	  
	    break;
		
	  case Moose::NeighborNeighbor:	
	  
	    if (vdotn >= 0.0) {
		 
          if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
		
            // negative gnd enters from neighbour		
		    jac += (-1.0) * ((vdotn * neigh_rho_gnd_neg * _u_neighbor[_qp]) 
			       / (_rho_tot_neighbor[_qp] * _rho_tot_neighbor[_qp])) 
		           * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
		
		  }
		
	    }

        if (vdotn < 0.0) {
		  
          if (_rho_tot_neighbor[_qp] > _rho_tot_tol) {
			
		    // positive gnd enters from neighbour
	        jac -= (-1.0) * ((vdotn * neigh_rho_gnd_pos * _u_neighbor[_qp]) 
			       / (_rho_tot_neighbor[_qp] * _rho_tot_neighbor[_qp])) 
	               * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];		
			
		  }		 

	    }	  	  
	  
	    break;
	}
	
  } else if (_rho_gnd_ot_coupled && jvar == _rho_gnd_ot_var) {
	  
    // TO DO
	jac += 0.0;
  
  } else {
	  
    jac += 0.0;  
	  
  }
  
  return jac; 
}
