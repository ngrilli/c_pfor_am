// Nicolò Grilli
// Università di Bristol
// 9 Novembre 2026

#include "CoupledCurlTestComponent.h"

registerMooseObject("c_pfor_amApp", CoupledCurlTestComponent);

InputParameters
CoupledCurlTestComponent::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Set the kernel variable to the weak form H * curl(test), where H is a vector "
                             "field made of three coupled variables Hx, Hy, Hz.");
  params.addRequiredCoupledVar("Hx", "Coupled variable representing the x-component of the vector field H");
  params.addRequiredCoupledVar("Hy", "Coupled variable representing the y-component of the vector field H");
  params.addRequiredCoupledVar("Hz", "Coupled variable representing the z-component of the vector field H");
  params.addRequiredParam<unsigned int>("component", "Component of the curl of Hx, Hy, Hz calculated by this kernel");
  return params;
}

CoupledCurlTestComponent::CoupledCurlTestComponent(const InputParameters & parameters)
  : Kernel(parameters),
    _Hx_var(coupled("Hx")),
    _Hy_var(coupled("Hy")),
    _Hz_var(coupled("Hz")),
    _Hx(coupledValue("Hx")),
    _Hy(coupledValue("Hy")),
    _Hz(coupledValue("Hz")),
    _component(getParam<unsigned int>("component"))
{
  if (_component >= LIBMESH_DIM)
    paramError("component", "Component too large for LIBMESH_DIM");
}

Real
CoupledCurlTestComponent::computeQpResidual()
{
  // Weak form: H * curl(test) = J * test
  // For component 0 (test function for Jx): test = (test_x, 0, 0)
  // curl(test) = (0, ∂test_x/∂z, -∂test_x/∂y)
  // H · curl(test) = Hy * ∂test_x/∂z - Hz * ∂test_x/∂y
  //
  // For component 1 (test function for Jy): test = (0, test_y, 0)
  // curl(test) = (-∂test_y/∂z, 0, ∂test_y/∂x)
  // H · curl(test) = -Hx * ∂test_y/∂z + Hz * ∂test_y/∂x
  //
  // For component 2 (test function for Jz): test = (0, 0, test_z)
  // curl(test) = (∂test_z/∂y, -∂test_z/∂x, 0)
  // H · curl(test) = Hx * ∂test_z/∂y - Hy * ∂test_z/∂x
  //
  // _grad_test[_i][_qp](0) = ∂test/∂x, _grad_test[_i][_qp](1) = ∂test/∂y, _grad_test[_i][_qp](2) = ∂test/∂z
  // where _i is the test function index for this kernel's variable, either Jx, Jy, or Jz, and _qp is the quadrature point index

  if (_component == 0) // test function for Jx component
    return _Hy[_qp] * _grad_test[_i][_qp](2) - _Hz[_qp] * _grad_test[_i][_qp](1);
  else if (_component == 1) // test function for Jy component
    return -_Hx[_qp] * _grad_test[_i][_qp](2) + _Hz[_qp] * _grad_test[_i][_qp](0);
  else if (_component == 2) // test function for Jz component
    return _Hx[_qp] * _grad_test[_i][_qp](1) - _Hy[_qp] * _grad_test[_i][_qp](0);
  else
    return 0.0;
}

Real
CoupledCurlTestComponent::computeQpJacobian()
{
  // The kernel variable u does not appear in the residual (H · curl(test) has no u dependence
  return 0.0;
}

Real
CoupledCurlTestComponent::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Derivatives with respect to the coupled variables Hx, Hy, Hz
  if (jvar == _Hx_var) {
    if (_component == 0)
      return 0.0; // Hx doesn't appear in component 0 residual
    else if (_component == 1)
      return -_phi[_j][_qp] * _grad_test[_i][_qp](2);
    else if (_component == 2)
      return _phi[_j][_qp] * _grad_test[_i][_qp](1);
  }
  else if (jvar == _Hy_var) {
    if (_component == 0)
      return _phi[_j][_qp] * _grad_test[_i][_qp](2);
    else if (_component == 1)
      return 0.0; // Hy doesn't appear in component 1 residual
    else if (_component == 2)
      return -_phi[_j][_qp] * _grad_test[_i][_qp](0);
  }
  else if (jvar == _Hz_var) {
    if (_component == 0)
      return -_phi[_j][_qp] * _grad_test[_i][_qp](1);
    else if (_component == 1)
      return _phi[_j][_qp] * _grad_test[_i][_qp](0);
    else if (_component == 2)
      return 0.0; // Hz doesn't appear in component 2 residual
  }
  return 0.0;
}