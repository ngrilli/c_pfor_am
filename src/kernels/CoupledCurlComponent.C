// Nicolò Grilli
// Università di Bristol
// 21 Ottobre 2026

#include "CoupledCurlComponent.h"

registerMooseObject("c_pfor_amApp", CoupledCurlComponent);

InputParameters
CoupledCurlComponent::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Set the kernel variable to a specified component of the curl of a vector made of " 
                             "three variables Hx, Hy, Hz, where all the three variables involved in the calculations "
                             "are considered coupled variables");
  params.addRequiredCoupledVar("Hx", "Coupled variable representing the x-component of the vector field");
  params.addRequiredCoupledVar("Hy", "Coupled variable representing the y-component of the vector field");
  params.addRequiredCoupledVar("Hz", "Coupled variable representing the z-component of the vector field");
  params.addRequiredParam<unsigned int>("component","Component of the curl of Hx, Hy, Hz calculated by this kernel");
  return params;
}

CoupledCurlComponent::CoupledCurlComponent(const InputParameters & parameters)
  : Kernel(parameters),
    _Hx_var(coupled("Hx")),
    _Hy_var(coupled("Hy")),
    _Hz_var(coupled("Hz")),
    _grad_Hx(coupledGradient("Hx")),
    _grad_Hy(coupledGradient("Hy")),
    _grad_Hz(coupledGradient("Hz")),
    _component(getParam<unsigned int>("component"))
{
  if (_component >= LIBMESH_DIM)
    paramError("component", "Component too large for LIBMESH_DIM");
}

Real
CoupledCurlComponent::computeQpResidual()
{
  if (_component == 0)
    return (_grad_Hz[_qp](1) - _grad_Hy[_qp](2)) * _test[_i][_qp];
  else if (_component == 1)
    return (_grad_Hx[_qp](2) - _grad_Hz[_qp](0)) * _test[_i][_qp];
  else if (_component == 2)
    return (_grad_Hy[_qp](0) - _grad_Hx[_qp](1)) * _test[_i][_qp];
  else
    return 0.0;
}

Real
CoupledCurlComponent::computeQpJacobian()
{
  return 0.0; // _u[_qp] is never used in the residual, so the Jacobian is always zero
}

Real
CoupledCurlComponent::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _Hx_var) {
    if (_component == 0)
      return 0.0;
    else if (_component == 1)
      return _grad_phi[_j][_qp](2) * _test[_i][_qp];
    else if (_component == 2)
      return -_grad_phi[_j][_qp](1) * _test[_i][_qp];
  }
  else if (jvar == _Hy_var) {
    if (_component == 0)
      return -_grad_phi[_j][_qp](2) * _test[_i][_qp];
    else if (_component == 1)
      return 0.0;
    else if (_component == 2)
      return _grad_phi[_j][_qp](0) * _test[_i][_qp];
  }
  else if (jvar == _Hz_var) {
    if (_component == 0)
      return _grad_phi[_j][_qp](1) * _test[_i][_qp];
    else if (_component == 1)
      return -_grad_phi[_j][_qp](0) * _test[_i][_qp];
    else if (_component == 2)
      return 0.0;
  }
  return 0.0;
}
