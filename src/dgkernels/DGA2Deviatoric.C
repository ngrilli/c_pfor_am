// Nicolo Grilli
// University of Bristol
// Daijun Hu
// National University of Singapore
// 22 Ottobre 2021

#include "DGA2Deviatoric.h"

registerMooseObject("MooseApp", DGA2Deviatoric);

InputParameters
DGA2Deviatoric::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription("Dislocation curvature diffusion "
                             "second term on the right hand side of equation (3) in "
							 "Stefan Sandfeld and Michael Zaiser (deviatoric part)"
							 "Pattern formation in a minimal model of continuum dislocation plasticity "
							 "Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp). "
							 "Discontinuous Galerkin formulation");
  params.addCoupledVar("rho_gnd_edge", 0.0, "Edge dislocation density: rho_x.");
  params.addCoupledVar("rho_gnd_screw", 0.0, "Screw dislocation density: rho_y.");
  params.addCoupledVar("dv_dx", 0.0, "Derivative of the velocity with respect to edge slip direction.");
  params.addCoupledVar("dv_dy", 0.0, "Derivative of the velocity with respect to screw slip direction.");
  params.addParam<Real>("dv_dx_max",1e9,"Max absolute value of dv_dx");
  params.addParam<Real>("dv_dy_max",1e9,"Max absolute value of dv_dy");
  params.addParam<Real>("ksabs_tol",0.000001,"Tolerance on small values of ksabs.");  
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");  
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  return params;
}

DGA2Deviatoric::DGA2Deviatoric(const InputParameters & parameters)
  : DGKernel(parameters),
    _rho_gnd_edge(coupledValue("rho_gnd_edge")), // Edge dislocation density: rho_x
    _rho_gnd_edge_coupled(isCoupled("rho_gnd_edge")),
    _rho_gnd_edge_var(_rho_gnd_edge_coupled ? coupled("rho_gnd_edge") : 0),
	_rho_gnd_edge_neighbor(coupledNeighborValue("rho_gnd_edge")),
    _rho_gnd_screw(coupledValue("rho_gnd_screw")), // Screw dislocation density: rho_y
    _rho_gnd_screw_coupled(isCoupled("rho_gnd_screw")),
    _rho_gnd_screw_var(_rho_gnd_screw_coupled ? coupled("rho_gnd_screw") : 0),
	_rho_gnd_screw_neighbor(coupledNeighborValue("rho_gnd_screw")),	
    _dv_dx(coupledValue("dv_dx")), // Derivative of the velocity with respect to edge slip direction
    _dv_dy(coupledValue("dv_dy")), // Derivative of the velocity with respect to screw slip direction
	_dv_dx_neighbor(coupledNeighborValue("dv_dx")), // same in the neighbouring element
	_dv_dy_neighbor(coupledNeighborValue("dv_dy")), // same in the neighbouring element
    _dv_dx_max(getParam<Real>("dv_dx_max")),
	_dv_dy_max(getParam<Real>("dv_dy_max")),
	_ksabs_tol(getParam<Real>("ksabs_tol")),
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
DGA2Deviatoric::getDislocationVelocity()
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
DGA2Deviatoric::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0.0;
  
  Real advected_quantity = 0.0; // term inside the d/dx or d/dy derivative
  Real neigh_advected_quantity = 0.0; // same in the neighbouring element
  
  Real ksabs; // modulus of k GND vector
  Real neigh_ksabs; // same in the neighbouring element
  
  Real diagterm; // diagonal term of deviatoric A2 matrix
  Real neigh_diagterm; // same in the neighbouring element
  
  Real outofdiagterm; // out of diagonal term of deviatoric A2 matrix  
  Real neigh_outofdiagterm; // same in the neighbouring element
  
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
  
  // sqrt(rho_x^2 + rho_y^2)
  ksabs = std::sqrt(_rho_gnd_edge[_qp]*_rho_gnd_edge[_qp]+_rho_gnd_screw[_qp]*_rho_gnd_screw[_qp]);
  neigh_ksabs = std::sqrt(_rho_gnd_edge_neighbor[_qp]*_rho_gnd_edge_neighbor[_qp]
              +_rho_gnd_screw_neighbor[_qp]*_rho_gnd_screw_neighbor[_qp]);
			  
  // rho_x^2 - rho_y^2
  diagterm = _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp] - _rho_gnd_screw[_qp]*_rho_gnd_screw[_qp];
  neigh_diagterm = _rho_gnd_edge_neighbor[_qp]*_rho_gnd_edge_neighbor[_qp]
                 - _rho_gnd_screw_neighbor[_qp]*_rho_gnd_screw_neighbor[_qp];
				 
  // rho_x rho_y
  outofdiagterm = _rho_gnd_edge[_qp]*_rho_gnd_screw[_qp];
  neigh_outofdiagterm = _rho_gnd_edge_neighbor[_qp]*_rho_gnd_screw_neighbor[_qp];
  
  if (ksabs > _ksabs_tol) {
    diagterm = 0.5 * diagterm / ksabs;
	outofdiagterm = outofdiagterm / ksabs;
  } else {
	diagterm = 0.0;  
	outofdiagterm = 0.0;
  }
  
  if (neigh_ksabs > _ksabs_tol) {
    neigh_diagterm = 0.5 * neigh_diagterm / neigh_ksabs;
	neigh_outofdiagterm = neigh_outofdiagterm / neigh_ksabs;
  } else {
	neigh_diagterm = 0.0;
	neigh_outofdiagterm = 0.0;
  }
  
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
	  advected_quantity -= diagterm * dv_dx;
	  advected_quantity -= outofdiagterm * dv_dy;
	  neigh_advected_quantity -= neigh_diagterm * neigh_dv_dx;
	  neigh_advected_quantity -= neigh_outofdiagterm * neigh_dv_dy;
	  break;
	case DisloCharacter::screw:
      advected_quantity += diagterm * dv_dy;
	  advected_quantity -= outofdiagterm * dv_dx;
	  neigh_advected_quantity += neigh_diagterm * neigh_dv_dy;
	  neigh_advected_quantity -= neigh_outofdiagterm * neigh_dv_dx;
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
DGA2Deviatoric::computeQpJacobian(Moose::DGJacobianType type)
{
  return 0.0;
}

// This term depend on the edge and screw dislocation densities
// So calculate derivative with respect to rho_x and rho_y
Real
DGA2Deviatoric::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real r = 0.0;
  
  Real advected_quantity = 0.0; // term inside the d/dx or d/dy derivative
  Real neigh_advected_quantity = 0.0; // same in the neighbouring element
  
  Real ksabs; // modulus of k GND vector
  Real neigh_ksabs; // same in the neighbouring element
  Real ksabs3; // modulus of k GND vector, power 3  
  Real neigh_ksabs3; // same in the neighbouring element

  // Derivatives with respect to _rho_gnd_edge and _rho_gnd_screw
  // of the diagonal term of deviatoric A2 matrix
  Real ddiagterm_drhoe;
  Real dneigh_diagterm_drhoe;
  Real ddiagterm_drhos; 
  Real dneigh_diagterm_drhos;

  // Derivatives with respect to _rho_gnd_edge and _rho_gnd_screw
  // of the out of diagonal term of deviatoric A2 matrix  
  Real doutofdiagterm_drhoe;
  Real dneigh_outofdiagterm_drhoe;
  Real doutofdiagterm_drhos;
  Real dneigh_outofdiagterm_drhos;
  
  // Velocity derivatives
  Real dv_dx;
  Real dv_dy;
  Real neigh_dv_dx;
  Real neigh_dv_dy;
  
  // Advected quantities for the Jacobian 
  Real d_advected = 0.0; // val
  Real d_neigh_advected = 0.0; // neigh_val;

  getDislocationVelocity();

  // scalar product between direction (edge or screw) 
  // and interface normal between this element and the neighbour
  Real vdotn;
  
  if (_rho_gnd_edge_coupled && jvar == _rho_gnd_edge_var)
  {	  
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
	  
    vdotn = _velocity * _normals[_qp];
	
    // sqrt(rho_x^2 + rho_y^2)
    ksabs = std::sqrt(_rho_gnd_edge[_qp]*_rho_gnd_edge[_qp]+_rho_gnd_screw[_qp]*_rho_gnd_screw[_qp]);
    neigh_ksabs = std::sqrt(_rho_gnd_edge_neighbor[_qp]*_rho_gnd_edge_neighbor[_qp]
                +_rho_gnd_screw_neighbor[_qp]*_rho_gnd_screw_neighbor[_qp]);
				
	// modulus of k GND vector, power 3
	// sqrt(rho_x^2 + rho_y^2)^3
    ksabs3 = std::pow(ksabs,3.0);
    neigh_ksabs3 = std::pow(neigh_ksabs,3.0);
	
    // Derivatives with respect to _rho_gnd_edge and _rho_gnd_screw
    // of the diagonal term of deviatoric A2 matrix
    ddiagterm_drhoe = _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp]
	                + 3.0*_rho_gnd_screw[_qp]*_rho_gnd_screw[_qp];
					
    ddiagterm_drhoe = 0.5 * _rho_gnd_edge[_qp] * ddiagterm_drhoe;
  
    dneigh_diagterm_drhoe = _rho_gnd_edge_neighbor[_qp]*_rho_gnd_edge_neighbor[_qp] 
	                      + 3.0*_rho_gnd_screw_neighbor[_qp]*_rho_gnd_screw_neighbor[_qp];
    dneigh_diagterm_drhoe = 0.5 * _rho_gnd_edge_neighbor[_qp] * dneigh_diagterm_drhoe;
	
	if (ksabs > _ksabs_tol) {
		
	  ddiagterm_drhoe = ddiagterm_drhoe / ksabs3;
      dneigh_diagterm_drhoe = dneigh_diagterm_drhoe / ksabs3;
	  
	} else {
		
	  ddiagterm_drhoe = 0.0;	
      dneigh_diagterm_drhoe = 0.0;
	  
	}	
	
    doutofdiagterm_drhoe = _rho_gnd_screw[_qp]*_rho_gnd_screw[_qp]*_rho_gnd_screw[_qp];
    dneigh_outofdiagterm_drhoe = _rho_gnd_screw_neighbor[_qp] * _rho_gnd_screw_neighbor[_qp] 
	                           * _rho_gnd_screw_neighbor[_qp];

	if (ksabs > _ksabs_tol) {
		
	  doutofdiagterm_drhoe = doutofdiagterm_drhoe / ksabs3;
      dneigh_outofdiagterm_drhoe = dneigh_outofdiagterm_drhoe / ksabs3;
	  
	} else {
		
	  doutofdiagterm_drhoe = 0.0;	
      dneigh_outofdiagterm_drhoe = 0.0;
	  
	}	
	
    switch (_dislo_character) // derivative type
    {
      case DisloCharacter::edge: // d/dx

        d_advected = (-1.0) * (ddiagterm_drhoe * dv_dx + doutofdiagterm_drhoe * dv_dy);
        d_neigh_advected = (-1.0) * (dneigh_diagterm_drhoe * neigh_dv_dx 
	                     + dneigh_outofdiagterm_drhoe * neigh_dv_dy);

        break;

      case DisloCharacter::screw: // d/dy

        d_advected = ddiagterm_drhoe * dv_dy - doutofdiagterm_drhoe * dv_dx;
        d_neigh_advected = dneigh_diagterm_drhoe * neigh_dv_dy 
	                   - dneigh_outofdiagterm_drhoe * neigh_dv_dx;

        break;

    } // end of switch case derivative type	
	
    switch (type) // Jacobian type
    {
      case Moose::ElementElement:
	  
        if (vdotn >= 0.0) {
          r += vdotn * d_advected *_phi[_j][_qp] * _test[_i][_qp];
        }
        break;

      case Moose::ElementNeighbor:

        if (vdotn < 0) {
          r += vdotn * d_neigh_advected * _phi_neighbor[_j][_qp] * _test[_i][_qp];
        }
        break;

      case Moose::NeighborElement:
      
        if (vdotn >= 0) {
          r -= vdotn * d_advected *_phi[_j][_qp] * _test_neighbor[_i][_qp];
        }
        break;
	    
      case Moose::NeighborNeighbor:
         
	    if (vdotn < 0) {
          r -= vdotn * d_neigh_advected *_phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
        }
        break;
		
    } // end of switch case Jacobian type	

  // end derivative with respect to rho_edge
	  
  } else if (_rho_gnd_screw_coupled && jvar == _rho_gnd_screw_var) {

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
	 
    vdotn = _velocity * _normals[_qp];
	
    // sqrt(rho_x^2 + rho_y^2)
    ksabs = std::sqrt(_rho_gnd_edge[_qp]*_rho_gnd_edge[_qp]+_rho_gnd_screw[_qp]*_rho_gnd_screw[_qp]);
    neigh_ksabs = std::sqrt(_rho_gnd_edge_neighbor[_qp]*_rho_gnd_edge_neighbor[_qp]
                +_rho_gnd_screw_neighbor[_qp]*_rho_gnd_screw_neighbor[_qp]);
				
	// modulus of k GND vector, power 3
	// sqrt(rho_x^2 + rho_y^2)^3
    ksabs3 = std::pow(ksabs,3.0);
    neigh_ksabs3 = std::pow(neigh_ksabs,3.0);
	
    // Derivatives with respect to _rho_gnd_edge and _rho_gnd_screw
    // of the diagonal term of deviatoric A2 matrix
    ddiagterm_drhos = _rho_gnd_screw[_qp]*_rho_gnd_screw[_qp]
	                + 3.0*_rho_gnd_edge[_qp]*_rho_gnd_edge[_qp];
					
    ddiagterm_drhos = (-0.5) * _rho_gnd_screw[_qp] * ddiagterm_drhos;
  
    dneigh_diagterm_drhos = _rho_gnd_screw_neighbor[_qp]*_rho_gnd_screw_neighbor[_qp] 
	                      + 3.0*_rho_gnd_edge_neighbor[_qp]*_rho_gnd_edge_neighbor[_qp];
    dneigh_diagterm_drhos = (-0.5) * _rho_gnd_screw_neighbor[_qp] * dneigh_diagterm_drhos;
	
	if (ksabs > _ksabs_tol) {
		
	  ddiagterm_drhos = ddiagterm_drhos / ksabs3;
      dneigh_diagterm_drhos = dneigh_diagterm_drhos / ksabs3;
	  
	} else {
		
	  ddiagterm_drhos = 0.0;	
      dneigh_diagterm_drhos = 0.0;
	  
	}

    doutofdiagterm_drhos = _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp]*_rho_gnd_edge[_qp];
    dneigh_outofdiagterm_drhos = _rho_gnd_edge_neighbor[_qp]*_rho_gnd_edge_neighbor[_qp]
	                           *_rho_gnd_edge_neighbor[_qp];

	if (ksabs > _ksabs_tol) {
		
	  doutofdiagterm_drhos = doutofdiagterm_drhos / ksabs3;
      dneigh_outofdiagterm_drhos = dneigh_outofdiagterm_drhos / ksabs3;
	  
	} else {
		
	  doutofdiagterm_drhos = 0.0;	
      dneigh_outofdiagterm_drhos = 0.0;
	
	}
	
    switch (_dislo_character) // derivative type
    {
      case DisloCharacter::edge: // d/dx

        d_advected = (-1.0) * (ddiagterm_drhos * dv_dx + doutofdiagterm_drhos * dv_dy);
        d_neigh_advected = (-1.0) * (dneigh_diagterm_drhos * neigh_dv_dx 
		                 + dneigh_outofdiagterm_drhos * neigh_dv_dy);

        break;

      case DisloCharacter::screw: // d/dy

        d_advected = ddiagterm_drhos * dv_dy - doutofdiagterm_drhos * dv_dx;
        d_neigh_advected = dneigh_diagterm_drhos * neigh_dv_dy 
		                 - dneigh_outofdiagterm_drhos * neigh_dv_dx;

        break;
		
    } // end of switch case derivative type	
	
    switch (type)
    {
      case Moose::ElementElement:
	  
        if (vdotn >= 0.0) {
          r += vdotn * d_advected *_phi[_j][_qp] * _test[_i][_qp];
        }
        break;

      case Moose::ElementNeighbor:

        if (vdotn < 0) {
          r += vdotn * d_neigh_advected * _phi_neighbor[_j][_qp] * _test[_i][_qp];
        }
        break;

      case Moose::NeighborElement:
      
        if (vdotn >= 0) {
          r -= vdotn * d_advected *_phi[_j][_qp] * _test_neighbor[_i][_qp];
        }
        break;
	    
      case Moose::NeighborNeighbor:
         
		if (vdotn < 0) {
          r -= vdotn * d_neigh_advected  *_phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
        }
        break;
	 
    } // end of switch case type	
  
  } // end derivative with respect to rho_screw

  return r;  
}
