// Nicolò Grilli
// University of Bristol
// 2 Luglio 2021

#include "DGAdvectionCoupledVConst.h"

registerMooseObject("MooseApp", DGAdvectionCoupledVConst);

InputParameters
DGAdvectionCoupledVConst::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription("DG upwinding for the advection of a coupled variable, constant V");
  params.addCoupledVar("rho_coupled", 0.0, "Coupled dislocation density in the flux term.");
  params.addRequiredParam<RealVectorValue>("velocity", "Velocity vector");
  return params;
}

DGAdvectionCoupledVConst::DGAdvectionCoupledVConst(const InputParameters & parameters)
  : DGKernel(parameters),
    _rho_coupled(coupledValue("rho_coupled")), // Coupled dislocation density in the flux term 
    _rho_coupled_coupled(isCoupled("rho_coupled")),
    _rho_coupled_var(_rho_coupled_coupled ? coupled("rho_coupled") : 0),
	_rho_neighbor(coupledNeighborValue("rho_coupled")),
	_velocity(getParam<RealVectorValue>("velocity"))
{
}

Real
DGAdvectionCoupledVConst::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  Real vdotn = _velocity * _normals[_qp];

  switch (type)
  {
    case Moose::Element:
      if (vdotn >= 0)
        r += vdotn * _rho_coupled[_qp] * _test[_i][_qp];
      else
        r += vdotn * _rho_neighbor[_qp] * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      if (vdotn >= 0)
        r -= vdotn * _rho_coupled[_qp] * _test_neighbor[_i][_qp];
      else
        r -= vdotn * _rho_neighbor[_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}

// Jacobian is zero because only coupled variable appears
Real DGAdvectionCoupledVConst::computeQpJacobian(Moose::DGJacobianType /*type*/) { return 0; }

Real
DGAdvectionCoupledVConst::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real r = 0;
  Real vdotn;

  if (_rho_coupled_coupled && jvar == _rho_coupled_var) {

    vdotn = _velocity * _normals[_qp];

    switch (type)
    {
      case Moose::ElementElement:
        if (vdotn >= 0)
          r += vdotn * _phi[_j][_qp] * _test[_i][_qp];
        break;

      case Moose::ElementNeighbor:
        if (vdotn < 0)
          r += vdotn * _phi_neighbor[_j][_qp] * _test[_i][_qp];
        break;

      case Moose::NeighborElement:
        if (vdotn >= 0)
          r -= vdotn * _phi[_j][_qp] * _test_neighbor[_i][_qp];
        break;

      case Moose::NeighborNeighbor:
        if (vdotn < 0)
          r -= vdotn * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
        break;
    }

  }

  return r;
}
