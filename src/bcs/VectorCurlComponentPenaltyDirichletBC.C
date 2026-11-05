// Nicolò Grilli
// Università di Bristol
// 5 Novembre 2026

#include "VectorCurlComponentPenaltyDirichletBC.h"
#include "Function.h"

registerMooseObject("c_pfor_amApp", VectorCurlComponentPenaltyDirichletBC);

InputParameters
VectorCurlComponentPenaltyDirichletBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<Real>("penalty", "The penalty coefficient");
  params.addRequiredParam<unsigned int>("component","Component of the curl of H to constrain (0=x, 1=y, 2=z)");
  params.addRequiredCoupledVar("H1", "Coupled variable whose derivative has positive sign in the curl component");
  params.addRequiredCoupledVar("H2", "Coupled variable whose derivative has negative sign in the curl component");
  params.addParam<FunctionName>("exact_curl_function", 0.0, "The function for the component of the exact curl");
  params.addClassDescription("Enforces a penalty Dirichlet boundary condition for a specified "
                             "component of the curl of a vector field H = (Hx, Hy, Hz), where H is "
                             "constructed from three coupled scalar variables. The BC penalizes the "
                             "difference between the computed curl component and a prescribed exact value.");
  return params;
}

VectorCurlComponentPenaltyDirichletBC::VectorCurlComponentPenaltyDirichletBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _penalty(getParam<Real>("penalty")),
    _component(getParam<unsigned int>("component")),
    _H1_var(coupled("H1")),
    _H2_var(coupled("H2")),
    _grad_H1(coupledGradient("H1")),
    _grad_H2(coupledGradient("H2")),
    _exact_curl_function(getFunction("exact_curl_function"))
{
  if (_component >= LIBMESH_DIM)
    paramError("component", "Component must be less than LIBMESH_DIM (0, 1, or 2)");
}

Real
VectorCurlComponentPenaltyDirichletBC::computeQpResidual()
{
  Real curl_H_component = computeCurlComponent();

  return _penalty * (curl_H_component - _exact_curl_function.value(_t, _q_point[_qp])) * _test[_i][_qp];
}

Real
VectorCurlComponentPenaltyDirichletBC::computeQpJacobian()
{
  // The residual does not depend on _u (the variable this BC is applied to), so the diagonal Jacobian is zero.
  return 0.0;
}

Real
VectorCurlComponentPenaltyDirichletBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0) // component 0: u=Hx, H1=Hz, H2=Hy -> _grad_H1[_qp](1) - _grad_H2[_qp](2)
  {
    if (jvar == _H1_var)
      return _grad_phi[_j][_qp](1) * _penalty * _test[_i][_qp];
    else if (jvar == _H2_var)
      return -_grad_phi[_j][_qp](2) * _penalty * _test[_i][_qp];
  }
  else if (_component == 1) // component 1: u=Hy, H1=Hx, H2=Hz -> _grad_H1[_qp](2) - _grad_H2[_qp](0)
  {
    if (jvar == _H1_var)
      return _grad_phi[_j][_qp](2) * _penalty * _test[_i][_qp];
    else if (jvar == _H2_var)
      return -_grad_phi[_j][_qp](0) * _penalty * _test[_i][_qp];
  }
  else if (_component == 2) // component 2: u=Hz, H1=Hy, H2=Hx -> _grad_H1[_qp](0) - _grad_H2[_qp](1)
  {
    if (jvar == _H1_var)
      return _grad_phi[_j][_qp](0) * _penalty * _test[_i][_qp];
    else if (jvar == _H2_var)
      return -_grad_phi[_j][_qp](1) * _penalty * _test[_i][_qp];
  }

  return 0.0;
}

Real
VectorCurlComponentPenaltyDirichletBC::computeCurlComponent() const
{
  // curl(H) = (dHz/dy - dHy/dz, dHx/dz - dHz/dx, dHy/dx - dHx/dy)
  // Cases:
  // component 0: u=Hx, H1=Hz, H2=Hy
  // component 1: u=Hy, H1=Hx, H2=Hz
  // component 2: u=Hz, H1=Hy, H2=Hx
  if (_component == 0)
    return _grad_H1[_qp](1) - _grad_H2[_qp](2);
  else if (_component == 1)
    return _grad_H1[_qp](2) - _grad_H2[_qp](0);
  else if (_component == 2)
    return _grad_H1[_qp](0) - _grad_H2[_qp](1);
  else
    return 0.0;
}