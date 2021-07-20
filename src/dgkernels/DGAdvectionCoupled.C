// Nicolò Grilli
// University of Bristol
// 30 Giugno 2021

#include "DGAdvectionCoupled.h"

registerMooseObject("MooseApp", DGAdvectionCoupled);

defineLegacyParams(DGAdvectionCoupled);

InputParameters
DGAdvectionCoupled::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription("DG upwinding for the advection of a coupled variable");
  params.addCoupledVar("rho_coupled", 0.0, "Coupled dislocation density in the flux term.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  params.addParam<bool>("is_edge_or_screw",false,"Check if is edge or screw, or total density");
  return params;
}

DGAdvectionCoupled::DGAdvectionCoupled(const InputParameters & parameters)
  : DGKernel(parameters),
    _rho_coupled(coupledValue("rho_coupled")), // Coupled dislocation density in the flux term 
    _rho_coupled_coupled(isCoupled("rho_coupled")),
    _rho_coupled_var(_rho_coupled_coupled ? coupled("rho_coupled") : 0),
	_rho_neighbor(coupledNeighborValue("rho_coupled")),
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
DGAdvectionCoupled::getDislocationVelocity()
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
DGAdvectionCoupled::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  
  getDislocationVelocity();

  Real vdotn = _velocity * _normals[_qp];
  Real rhocvdotn = vdotn * _rho_coupled[_qp];
  Real rhovdotn = vdotn * _u[_qp];
  Real neighrhocvdotn = vdotn * _rho_neighbor[_qp];

  switch (type)
  {
    case Moose::Element:	
	  
	  if (_is_edge_or_screw) { // Case with rho_edge or rho_screw = _u[_qp]
	  
	    if (rhovdotn >= 0)
          r += rhocvdotn * _test[_i][_qp];
        else
          r += neighrhocvdotn * _test[_i][_qp];
		  
	  } else { // Case with rho_edge or rho_screw = _rho_coupled[_qp]
		  
	    if (rhocvdotn >= 0)
          r += rhocvdotn * _test[_i][_qp];
        else
          r += neighrhocvdotn * _test[_i][_qp];	
	  
	  }

      break;

    case Moose::Neighbor:
	
	  if (_is_edge_or_screw) {
		  
	    if (rhovdotn >= 0)
          r -= rhocvdotn * _test_neighbor[_i][_qp];
        else
          r -= neighrhocvdotn * _test_neighbor[_i][_qp];
		  
	  } else {
		  
        if (rhocvdotn >= 0)
          r -= rhocvdotn * _test_neighbor[_i][_qp];
        else
          r -= neighrhocvdotn * _test_neighbor[_i][_qp];	
	  
	  }
	  
      break;
  }

  return r;
}

// Jacobian is zero because only coupled variable appears
Real DGAdvectionCoupled::computeQpJacobian(Moose::DGJacobianType /*type*/) { return 0; }

Real
DGAdvectionCoupled::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real r = 0;
  Real vdotn;
  Real rhocvdotn;
  Real rhovdotn;

  if (_rho_coupled_coupled && jvar == _rho_coupled_var) {
	  
    getDislocationVelocity();

    vdotn = _velocity * _normals[_qp];
	rhocvdotn = vdotn * _rho_coupled[_qp];
	rhovdotn = vdotn * _u[_qp];

    switch (type)
    {
      case Moose::ElementElement:
	  
	    if (_is_edge_or_screw) {
			
	      if (rhovdotn >= 0)
            r += vdotn * _phi[_j][_qp] * _test[_i][_qp];
			
		} else {
			
		  if (rhocvdotn >= 0)
            r += vdotn * _phi[_j][_qp] * _test[_i][_qp];
	  
		}
		
        break;

      case Moose::ElementNeighbor:
	  
	    if (_is_edge_or_screw) {
			
	      if (rhovdotn < 0)
            r += vdotn * _phi_neighbor[_j][_qp] * _test[_i][_qp];
			
		} else {
			
          if (rhocvdotn < 0)
            r += vdotn * _phi_neighbor[_j][_qp] * _test[_i][_qp];	
		
		}
		
        break;

      case Moose::NeighborElement:
	  
	    if (_is_edge_or_screw) {
			
	      if (rhovdotn >= 0)
            r -= vdotn * _phi[_j][_qp] * _test_neighbor[_i][_qp];
			
		} else {
			
          if (rhocvdotn >= 0)
            r -= vdotn * _phi[_j][_qp] * _test_neighbor[_i][_qp];
		
		}
		
        break;

      case Moose::NeighborNeighbor:
	  
        if (_is_edge_or_screw) {
			
	      if (rhovdotn < 0)
            r -= vdotn * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
			
		} else {
			
          if (rhocvdotn < 0)
            r -= vdotn * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];	
		
		}
		
        break;

    }

  }

  return r;
}
