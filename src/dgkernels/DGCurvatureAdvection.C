// Nicolo Grilli
// University of Bristol
// Daijun Hu
// National University of Singapore
// 26 Luglio 2021

#include "DGCurvatureAdvection.h"

registerMooseObject("MooseApp", DGCurvatureAdvection);

InputParameters
DGCurvatureAdvection::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription("DG Dislocation curvature advection."
							 "Upwind condition is calculated both on edge/screw dislocations "
							 "in this element and on the neighbouring element.");
  params.addCoupledVar("rho_gnd", 0.0, "GND dislocation density: rho_x or rho_y for edge or screw.");
  params.addCoupledVar("rho_tot", 0.0, "Total dislocation density: rho_t.");
  params.addParam<Real>("rho_tot_tol",0.000001,"Tolerance on small values of rho_tot.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  return params;
}

DGCurvatureAdvection::DGCurvatureAdvection(const InputParameters & parameters)
  : DGKernel(parameters),
    _rho_gnd(coupledValue("rho_gnd")), // GND dislocation density: rho_x or rho_y for edge or screw
    _rho_gnd_coupled(isCoupled("rho_gnd")),
    _rho_gnd_var(_rho_gnd_coupled ? coupled("rho_gnd") : 0),
    _rho_tot(coupledValue("rho_tot")), // Total dislocation density: rho_t
    _rho_tot_coupled(isCoupled("rho_tot")),
    _rho_tot_var(_rho_tot_coupled ? coupled("rho_tot") : 0),
	_rho_gnd_neighbor(coupledNeighborValue("rho_gnd")),
    _rho_tot_neighbor(coupledNeighborValue("rho_tot")),
	_rho_tot_tol(getParam<Real>("rho_tot_tol")),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
    _dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")), // Velocity value (signed)
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>())
{
}

// read dislocation velocity from material object
// and store in _velocity
void
DGCurvatureAdvection::getDislocationVelocity()
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
DGCurvatureAdvection::computeQpResidual(Moose::DGResidualType type)
{ 
  Real r = 0; // Residual for output
  
  getDislocationVelocity();

  Real vdotn = _velocity * _normals[_qp];
  
  Real rho_gnd_vdotn = vdotn * _rho_gnd[_qp];
  Real neigh_rho_gnd_vdotn = vdotn * _rho_gnd_neighbor[_qp];

  switch (type)
  {
    case Moose::Element:
	  	  
	  if (rho_gnd_vdotn >= 0.0 && _rho_tot[_qp] > _rho_tot_tol)
        r += ((_rho_gnd[_qp] * _u[_qp]) / _rho_tot[_qp]) 
	         * _test[_i][_qp];
	    
      if (neigh_rho_gnd_vdotn < 0.0 && _rho_tot_neighbor[_qp] > _rho_tot_tol)
	    r += ((_rho_gnd_neighbor[_qp] * _u_neighbor[_qp]) / _rho_tot_neighbor[_qp]) 
	         * _test[_i][_qp];	
		  
      break;

    case Moose::Neighbor:
	  		  
	  if (rho_gnd_vdotn >= 0.0 && _rho_tot[_qp] > _rho_tot_tol)
        r -= ((_rho_gnd[_qp] * _u[_qp]) / _rho_tot[_qp]) 
	         * _test_neighbor[_i][_qp];
	  
	  if (neigh_rho_gnd_vdotn < 0.0 && _rho_tot_neighbor[_qp] > _rho_tot_tol)
	    r -= ((_rho_gnd_neighbor[_qp] * _u_neighbor[_qp]) / _rho_tot_neighbor[_qp])  
	         * _test_neighbor[_i][_qp];
		  
      break;
  }

  return r;
}


Real 
DGCurvatureAdvection::computeQpJacobian(Moose::DGJacobianType type) 
{ 
  Real jac = 0; // Jacobian for output
  
  getDislocationVelocity();
  
  Real vdotn = _velocity * _normals[_qp];
  
  Real rho_gnd_vdotn = vdotn * _rho_gnd[_qp];
  Real neigh_rho_gnd_vdotn = vdotn * _rho_gnd_neighbor[_qp];

  switch (type)
  {
    case Moose::ElementElement:
	  	  
	  if (rho_gnd_vdotn >= 0.0 && _rho_tot[_qp] > _rho_tot_tol)
        jac += ((_rho_gnd[_qp] * _phi[_j][_qp]) / _rho_tot[_qp]) 
	           * _test[_i][_qp];	
		  
      break;

    case Moose::ElementNeighbor:
	    
      if (neigh_rho_gnd_vdotn < 0.0 && _rho_tot_neighbor[_qp] > _rho_tot_tol)
	    jac += ((_rho_gnd_neighbor[_qp] * _phi_neighbor[_j][_qp]) / _rho_tot_neighbor[_qp]) 
	           * _test[_i][_qp];	
		  
      break;
	  
	case Moose::NeighborElement:
	
	  if (rho_gnd_vdotn >= 0.0 && _rho_tot[_qp] > _rho_tot_tol)
        jac -= ((_rho_gnd[_qp] * _phi[_j][_qp]) / _rho_tot[_qp]) 
	           * _test_neighbor[_i][_qp];
	
	  break;
	
	case Moose::NeighborNeighbor:
	
	  if (neigh_rho_gnd_vdotn < 0.0 && _rho_tot_neighbor[_qp] > _rho_tot_tol)
	    jac -= ((_rho_gnd_neighbor[_qp] * _phi_neighbor[_j][_qp]) / _rho_tot_neighbor[_qp])  
	           * _test_neighbor[_i][_qp];	
	
	  break;
  }  

  return jac; 
}

Real
DGCurvatureAdvection::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real jac = 0; // Jacobian for output
  
  getDislocationVelocity();
  
  Real vdotn = _velocity * _normals[_qp];
  
  Real rho_gnd_vdotn = vdotn * _rho_gnd[_qp];
  Real neigh_rho_gnd_vdotn = vdotn * _rho_gnd_neighbor[_qp];
  
  if (_rho_gnd_coupled && jvar == _rho_gnd_var) {
	  
    // derivative with respect to _rho_gnd and _rho_gnd_neighbor
	
	switch (type)
    {
      case Moose::ElementElement:
	  
	    if (rho_gnd_vdotn >= 0.0 && _rho_tot[_qp] > _rho_tot_tol)
          jac += ((_phi[_j][_qp] * _u[_qp]) / _rho_tot[_qp]) 
	             * _test[_i][_qp];	  
		
        break;

      case Moose::ElementNeighbor:
	  
        if (neigh_rho_gnd_vdotn < 0.0 && _rho_tot_neighbor[_qp] > _rho_tot_tol)
	      jac += ((_phi_neighbor[_j][_qp] * _u_neighbor[_qp]) / _rho_tot_neighbor[_qp]) 
	             * _test[_i][_qp];	  
		
        break;	  
	  
	  case Moose::NeighborElement:
	  
	    if (rho_gnd_vdotn >= 0.0 && _rho_tot[_qp] > _rho_tot_tol)
          jac -= ((_phi[_j][_qp] * _u[_qp]) / _rho_tot[_qp]) 
	           * _test_neighbor[_i][_qp];	  
		
        break;

	  case Moose::NeighborNeighbor:
	  
	    if (neigh_rho_gnd_vdotn < 0.0 && _rho_tot_neighbor[_qp] > _rho_tot_tol)
	      jac -= ((_phi_neighbor[_j][_qp] * _u_neighbor[_qp]) / _rho_tot_neighbor[_qp])  
	           * _test_neighbor[_i][_qp];	  
		
        break;	  
	}
  
  } else if (_rho_tot_coupled && jvar == _rho_tot_var) {
	  
    // derivative with respect to _rho_tot and _rho_tot_neighbor
	switch (type)
    {
	  case Moose::ElementElement:
	  
	    if (rho_gnd_vdotn >= 0.0 && _rho_tot[_qp] > _rho_tot_tol)
          jac += (-1.0) * ((_rho_gnd[_qp] * _u[_qp]) / (_rho_tot[_qp] * _rho_tot[_qp])) 
	             * _phi[_j][_qp] * _test[_i][_qp];	  
	  
	    break;
		
	  case Moose::ElementNeighbor:
	  
        if (neigh_rho_gnd_vdotn < 0.0 && _rho_tot_neighbor[_qp] > _rho_tot_tol)
	      jac += (-1.0) * ((_rho_gnd_neighbor[_qp] * _u_neighbor[_qp]) 
	             / (_rho_tot_neighbor[_qp] * _rho_tot_neighbor[_qp])) 
	             * _phi_neighbor[_j][_qp] * _test[_i][_qp];		  
	  
	    break;
		
	  case Moose::NeighborElement:
	  
	    if (rho_gnd_vdotn >= 0.0 && _rho_tot[_qp] > _rho_tot_tol)
          jac -= (-1.0) * ((_rho_gnd[_qp] * _u[_qp]) / (_rho_tot[_qp] * _rho_tot[_qp])) 
	             * _phi[_j][_qp] * _test_neighbor[_i][_qp];	  
	  
	    break;
		
	  case Moose::NeighborNeighbor:	
	  
	    if (neigh_rho_gnd_vdotn < 0.0 && _rho_tot_neighbor[_qp] > _rho_tot_tol)
	      jac -= (-1.0) * ((_rho_gnd_neighbor[_qp] * _u_neighbor[_qp]) 
	             / (_rho_tot_neighbor[_qp] * _rho_tot_neighbor[_qp]))  
	             * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];	  
	  
	    break;
	}
	
  } else {
	  
    jac = 0.0;  
	  
  }
  
  return jac; 
}
