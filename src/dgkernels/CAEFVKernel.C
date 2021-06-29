// Nicolò Grilli
// University of Bristol
// 29 Giugno 2021

#include "CAEFVKernel.h"

registerMooseObject("RdgApp", CAEFVKernel);

InputParameters
CAEFVKernel::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription(
      "A dgkernel for the coupled advection equation using a cell-centered finite volume method.");
  MooseEnum component("concentration");
  params.addParam<MooseEnum>("component", component, "Choose one of the equations");
  params.addRequiredCoupledVar("rho_coupled", "Coupled dislocation density in the flux term");
  params.addRequiredParam<UserObjectName>("flux", "Name of the internal side flux object to use");
  return params;
}

CAEFVKernel::CAEFVKernel(const InputParameters & parameters)
  : DGKernel(parameters),
    _component(getParam<MooseEnum>("component")),
    _rho_coupled_c1(coupledValue("rho_coupled")),
    _rho_coupled_c2(coupledNeighborValue("rho_coupled")),
    _rho_coupled_1(getMaterialProperty<Real>("rho_coupled")),
    _rho_coupled_2(getNeighborMaterialProperty<Real>("rho_coupled")),
    _flux(getUserObject<InternalSideFluxBase>("flux"))
{
}

CAEFVKernel::~CAEFVKernel() {}

Real
CAEFVKernel::computeQpResidual(Moose::DGResidualType type)
{
  // assemble the input vectors, which are
  //   the reconstructed linear monomial
  //   extrapolated at side center from the current and neighbor elements
  std::vector<Real> uvec1 = {_rho_coupled_1[_qp]};
  std::vector<Real> uvec2 = {_rho_coupled_2[_qp]};

  // calculate the flux
  const auto & flux = _flux.getFlux(
      _current_side, _current_elem->id(), _neighbor_elem->id(), uvec1, uvec2, _normals[_qp]);

  // distribute the contribution to the current and neighbor elements
  switch (type)
  {
    case Moose::Element:
      return flux[_component] * _test[_i][_qp];

    case Moose::Neighbor:
      return -flux[_component] * _test_neighbor[_i][_qp];
  }

  return 0.0;
}

Real
CAEFVKernel::computeQpJacobian(Moose::DGJacobianType type)
{
  // assemble the input vectors, which are
  //   the constant monomial from the current and neighbor elements
  std::vector<Real> uvec1 = {_rho_coupled_c1[_qp]};
  std::vector<Real> uvec2 = {_rho_coupled_c2[_qp]};

  // calculate the Jacobian matrices
  const auto & fjac1 = _flux.getJacobian(Moose::Element,
                                         _current_side,
                                         _current_elem->id(),
                                         _neighbor_elem->id(),
                                         uvec1,
                                         uvec2,
                                         _normals[_qp]);

  const auto & fjac2 = _flux.getJacobian(Moose::Neighbor,
                                         _current_side,
                                         _current_elem->id(),
                                         _neighbor_elem->id(),
                                         uvec1,
                                         uvec2,
                                         _normals[_qp]);

  // distribute the contribution to the current and neighbor elements
  switch (type)
  {
    case Moose::ElementElement:
      return fjac1(_component, _component) * _phi[_j][_qp] * _test[_i][_qp];

    case Moose::ElementNeighbor:
      return fjac2(_component, _component) * _phi_neighbor[_j][_qp] * _test[_i][_qp];

    case Moose::NeighborElement:
      return -fjac1(_component, _component) * _phi[_j][_qp] * _test_neighbor[_i][_qp];

    case Moose::NeighborNeighbor:
      return -fjac2(_component, _component) * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
  }

  return 0.0;
}
