// Nicolò Grilli
// Università di Bristol
// 12 Settembre 2026

#include "CurlComponent.h"

registerMooseObject("c_pfor_amApp", CurlComponent);

InputParameters
CurlComponent::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Set the kernel variable to a specified component of the curl of a vector made of " 
                             "three variables Hx, Hy, Hz, where the two variables involved in the calculations "
                             "are considered coupled variables, while the variable not involved is considered "
                             "an independent variable");
  params.addRequiredCoupledVar("H1", "Coupled variable whose derivative has positive sign in the curl component");
  params.addRequiredCoupledVar("H2", "Coupled variable whose derivative has negative sign in the curl component");
  params.addRequiredParam<unsigned int>("component","Component of the curl of Hx, Hy, Hz calculated by this kernel");
  return params;
}

CurlComponent::CurlComponent(const InputParameters & parameters)
  : Kernel(parameters),
    _H1_var(coupled("H1")),
    _H2_var(coupled("H2")),
    _grad_H1(coupledGradient("H1")),
    _grad_H2(coupledGradient("H2")),
    _component(getParam<unsigned int>("component"))
{
  if (_component >= LIBMESH_DIM)
    paramError("component", "Component too large for LIBMESH_DIM");
}

Real
CurlComponent::computeQpResidual()
{
  if (_component == 0) // H1 = Hz, H2 = Hy
    return (_grad_H1[_qp](1) - _grad_H2[_qp](2)) * _test[_i][_qp];
  else if (_component == 1) // H1 = Hx, H2 = Hz
    return (_grad_H1[_qp](2) - _grad_H2[_qp](0)) * _test[_i][_qp];
  else if (_component == 2) // H1 = Hy, H2 = Hx
    return (_grad_H1[_qp](0) - _grad_H2[_qp](1)) * _test[_i][_qp];
  else
    return 0.0;
}

Real
CurlComponent::computeQpJacobian()
{
  return 0.0; // _u[_qp] is never used in the residual, so the Jacobian is always zero
}

Real
CurlComponent::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _H1_var) {
    if (_component == 0) // H1 = Hz, H2 = Hy
      return _grad_phi[_j][_qp](1) * _test[_i][_qp];
    else if (_component == 1) // H1 = Hx, H2 = Hz
      return _grad_phi[_j][_qp](2) * _test[_i][_qp];
    else if (_component == 2) // H1 = Hy, H2 = Hx
      return _grad_phi[_j][_qp](0) * _test[_i][_qp];
  }
  else if (jvar == _H2_var) {
    if (_component == 0) // H1 = Hz, H2 = Hy
      return -_grad_phi[_j][_qp](2) * _test[_i][_qp];
    else if (_component == 1) // H1 = Hx, H2 = Hz
      return -_grad_phi[_j][_qp](0) * _test[_i][_qp];
    else if (_component == 2) // H1 = Hy, H2 = Hx
      return -_grad_phi[_j][_qp](1) * _test[_i][_qp];
  } 
  return 0.0;
}
