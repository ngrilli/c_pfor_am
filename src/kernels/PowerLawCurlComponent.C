// Nicolò Grilli
// Università di Bristol
// 13 Settembre 2026

#include "PowerLawCurlComponent.h"

registerMooseObject("c_pfor_amApp", PowerLawCurlComponent);

InputParameters
PowerLawCurlComponent::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Set the kernel variable to a specified component of the curl of a vector made of " 
                             "three variables Jx, Jy, Jz, multiplied by a prefactor power law of the magnitude of the vector (Jx,Jy,Jz). "
                             "Two of the variables involved in the calculations are considered coupled variables, "
                             "while the third variable is considered an independent variable");
  params.addRequiredCoupledVar("J1", "Coupled variable whose derivative has positive sign in the curl component");
  params.addRequiredCoupledVar("J2", "Coupled variable whose derivative has negative sign in the curl component");
  params.addRequiredParam<unsigned int>("component","Component of the curl of (Jx,Jy,Jz), multiplied by the power law, calculated by this kernel");
  params.addRequiredParam<Real>("E0","Critical electric field");
  params.addRequiredParam<Real>("J0","Critical current density");
  params.addRequiredParam<Real>("n","Power-law exponent");
  return params;
}

PowerLawCurlComponent::PowerLawCurlComponent(const InputParameters & parameters)
  : Kernel(parameters),
    _J1_var(coupled("J1")),
    _J2_var(coupled("J2")),
    _J1(coupledValue("J1")),
    _J2(coupledValue("J2")),
    _grad_J1(coupledGradient("J1")),
    _grad_J2(coupledGradient("J2")),
    _component(getParam<unsigned int>("component")),
    _E0(getParam<Real>("E0")),
    _J0(getParam<Real>("J0")),
    _n(getParam<Real>("n"))
{
  if (_component >= LIBMESH_DIM)
    paramError("component", "Component too large for LIBMESH_DIM");
}

Real
PowerLawCurlComponent::computeQpResidual()
{
  const Real Jmag = std::sqrt(_J1[_qp] * _J1[_qp] + _J2[_qp] * _J2[_qp] + _u[_qp] * _u[_qp]);

  // Avoid division by zero in the power-law expression
  if (Jmag == 0.0)
    return 0.0;

  // Power law prefactor 
  const Real f = (_E0 / _J0) * std::pow(Jmag / _J0, _n - 1.0);

  // Spatial derivatives of |J| = Jmag
  const Real dJmag_dx = (_J1[_qp] * _grad_J1[_qp](0) + _J2[_qp] * _grad_J2[_qp](0) + _u[_qp] * _grad_u[_qp](0)) / Jmag;
  const Real dJmag_dy = (_J1[_qp] * _grad_J1[_qp](1) + _J2[_qp] * _grad_J2[_qp](1) + _u[_qp] * _grad_u[_qp](1)) / Jmag;
  const Real dJmag_dz = (_J1[_qp] * _grad_J1[_qp](2) + _J2[_qp] * _grad_J2[_qp](2) + _u[_qp] * _grad_u[_qp](2)) / Jmag;

  // Spatial derivatives of prefactor f
  const Real df_dJmag = f * (_n - 1.0) / Jmag;
  const Real df_dx = df_dJmag * dJmag_dx;
  const Real df_dy = df_dJmag * dJmag_dy;
  const Real df_dz = df_dJmag * dJmag_dz;

  if (_component == 0) // J1 = Jz, J2 = Jy, u = Jx
    return (f * (_grad_J1[_qp](1) - _grad_J2[_qp](2)) + _J1[_qp] * df_dy - _J2[_qp] * df_dz) * _test[_i][_qp];
  else if (_component == 1) // J1 = Jx, J2 = Jz, u = Jy
    return (f * (_grad_J1[_qp](2) - _grad_J2[_qp](0)) + _J1[_qp] * df_dz - _J2[_qp] * df_dx) * _test[_i][_qp];
  else if (_component == 2) // J1 = Jy, J2 = Jx, u = Jz
    return (f * (_grad_J1[_qp](0) - _grad_J2[_qp](1)) + _J1[_qp] * df_dx - _J2[_qp] * df_dy) * _test[_i][_qp];
  else
    return 0.0;
}

Real
PowerLawCurlComponent::computeQpJacobian()
{
  const Real Jmag = std::sqrt(_J1[_qp] * _J1[_qp] + _J2[_qp] * _J2[_qp] + _u[_qp] * _u[_qp]);

  // Avoid division by zero in the power-law expression
  if (Jmag == 0.0)
    return 0.0;

  // Power law prefactor and its derivative with respect to Jmag
  const Real f = (_E0 / _J0) * std::pow(Jmag / _J0, _n - 1.0);
  const Real df_dJmag = f * (_n - 1.0) / Jmag;

  // Derivative of f with respect to variable _u
  const Real dJmag_du = (_u[_qp] * _phi[_j][_qp]) / Jmag;
  const Real df_du = df_dJmag * dJmag_du * _phi[_j][_qp];

  // Spatial derivatives of |J| = Jmag
  const Real dJmag_numerator_dx = _J1[_qp] * _grad_J1[_qp](0) + _J2[_qp] * _grad_J2[_qp](0) + _u[_qp] * _grad_u[_qp](0);
  const Real dJmag_numerator_dy = _J1[_qp] * _grad_J1[_qp](1) + _J2[_qp] * _grad_J2[_qp](1) + _u[_qp] * _grad_u[_qp](1);
  const Real dJmag_numerator_dz = _J1[_qp] * _grad_J1[_qp](2) + _J2[_qp] * _grad_J2[_qp](2) + _u[_qp] * _grad_u[_qp](2);
  const Real dJmag_dx = dJmag_numerator_dx / Jmag;
  const Real dJmag_dy = dJmag_numerator_dy / Jmag;
  const Real dJmag_dz = dJmag_numerator_dz / Jmag;

  // Spatial derivatives of f
  const Real df_dx = df_dJmag * dJmag_dx;
  const Real df_dy = df_dJmag * dJmag_dy;
  const Real df_dz = df_dJmag * dJmag_dz;

  // Derivative of df_dJmag, df_dx, df_dy, df_dz with respect to u
  const Real d2f_dJmag_du = (_n - 1.0) * (df_du / Jmag - (f * dJmag_du) / (Jmag * Jmag)) * _phi[_j][_qp];
  const Real d2f_dx_du = d2f_dJmag_du * dJmag_dx + df_dJmag * ((_grad_u[_qp](0) * _phi[_j][_qp] + _u[_qp] * _grad_phi[_j][_qp](0)) / Jmag - (dJmag_numerator_dx * dJmag_du) / (Jmag * Jmag));
  const Real d2f_dy_du = d2f_dJmag_du * dJmag_dy + df_dJmag * ((_grad_u[_qp](1) * _phi[_j][_qp] + _u[_qp] * _grad_phi[_j][_qp](1)) / Jmag - (dJmag_numerator_dy * dJmag_du) / (Jmag * Jmag));
  const Real d2f_dz_du = d2f_dJmag_du * dJmag_dz + df_dJmag * ((_grad_u[_qp](2) * _phi[_j][_qp] + _u[_qp] * _grad_phi[_j][_qp](2)) / Jmag - (dJmag_numerator_dz * dJmag_du) / (Jmag * Jmag));

  if (_component == 0) // J1 = Jz, J2 = Jy, u = Jx
    return (df_du * (_grad_J1[_qp](1) - _grad_J2[_qp](2)) + _J1[_qp] * d2f_dy_du - _J2[_qp] * d2f_dz_du) * _test[_i][_qp];
  else if (_component == 1) // J1 = Jx, J2 = Jz, u = Jy
    return (df_du * (_grad_J1[_qp](2) - _grad_J2[_qp](0)) + _J1[_qp] * d2f_dz_du - _J2[_qp] * d2f_dx_du) * _test[_i][_qp];
  else if (_component == 2) // J1 = Jy, J2 = Jx, u = Jz
    return (df_du * (_grad_J1[_qp](0) - _grad_J2[_qp](1)) + _J1[_qp] * d2f_dx_du - _J2[_qp] * d2f_dy_du) * _test[_i][_qp];
  else
    return 0.0;
}

Real
PowerLawCurlComponent::computeQpOffDiagJacobian(unsigned int jvar) // TO MODIFY
{
  if (jvar == _J1_var) {
    if (_component == 0) // H1 = Hz, H2 = Hy
      return _grad_phi[_j][_qp](1) * _test[_i][_qp];
    else if (_component == 1) // H1 = Hx, H2 = Hz
      return _grad_phi[_j][_qp](2) * _test[_i][_qp];
    else if (_component == 2) // H1 = Hy, H2 = Hx
      return _grad_phi[_j][_qp](0) * _test[_i][_qp];
  }
  else if (jvar == _J2_var) {
    if (_component == 0) // H1 = Hz, H2 = Hy
      return -_grad_phi[_j][_qp](2) * _test[_i][_qp];
    else if (_component == 1) // H1 = Hx, H2 = Hz
      return -_grad_phi[_j][_qp](0) * _test[_i][_qp];
    else if (_component == 2) // H1 = Hy, H2 = Hx
      return -_grad_phi[_j][_qp](1) * _test[_i][_qp];
  } 
  return 0.0;
}
