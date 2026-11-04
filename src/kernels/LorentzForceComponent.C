// Nicolò Grilli
// Università di Bristol
// 4 Novembre 2026

#include "LorentzForceComponent.h"

registerMooseObject("c_pfor_amApp", LorentzForceComponent);

InputParameters
LorentzForceComponent::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates one component of the magnetic Lorentz force density: "
                             "f = mu_0 * (J x H), where J is the current density vector and "
                             "H is the magnetic field strength vector. The force is per unit volume.");
  params.addRequiredCoupledVar("Hx", "Coupled variable for x-component of magnetic field strength");
  params.addRequiredCoupledVar("Hy", "Coupled variable for y-component of magnetic field strength");
  params.addRequiredCoupledVar("Hz", "Coupled variable for z-component of magnetic field strength");
  params.addRequiredCoupledVar("Jx", "Coupled variable for x-component of current density");
  params.addRequiredCoupledVar("Jy", "Coupled variable for y-component of current density");
  params.addRequiredCoupledVar("Jz", "Coupled variable for z-component of current density");
  params.addRequiredParam<unsigned int>("component", "Component of the Lorentz force vector to compute (0=x, 1=y, 2=z)");
  params.addParam<Real>("mu_0", 1.25663706212e-6, "Magnetic permeability");
  return params;
}

LorentzForceComponent::LorentzForceComponent(const InputParameters & parameters)
  : Kernel(parameters),
    _Hx_var(coupled("Hx")),
    _Hy_var(coupled("Hy")),
    _Hz_var(coupled("Hz")),
    _Jx_var(coupled("Jx")),
    _Jy_var(coupled("Jy")),
    _Jz_var(coupled("Jz")),
    _Hx(coupledValue("Hx")),
    _Hy(coupledValue("Hy")),
    _Hz(coupledValue("Hz")),
    _Jx(coupledValue("Jx")),
    _Jy(coupledValue("Jy")),
    _Jz(coupledValue("Jz")),
    _component(getParam<unsigned int>("component")),
    _mu_0(getParam<Real>("mu_0"))
{
  if (_component >= LIBMESH_DIM)
    paramError("component", "Component too large for LIBMESH_DIM");
}

Real
LorentzForceComponent::computeQpResidual()
{
  // Lorentz force density: f = mu_0 * (J x H)
  // Component 0 (x): f_x = mu_0 * (Jy * Hz - Jz * Hy)
  // Component 1 (y): f_y = mu_0 * (Jz * Hx - Jx * Hz)
  // Component 2 (z): f_z = mu_0 * (Jx * Hy - Jy * Hx)

  Real force = 0.0;

  if (_component == 0) // x-component
    force = _Jy[_qp] * _Hz[_qp] - _Jz[_qp] * _Hy[_qp];
  else if (_component == 1) // y-component
    force = _Jz[_qp] * _Hx[_qp] - _Jx[_qp] * _Hz[_qp];
  else if (_component == 2) // z-component
    force = _Jx[_qp] * _Hy[_qp] - _Jy[_qp] * _Hx[_qp];

  return _mu_0 * force * _test[_i][_qp];
}

Real
LorentzForceComponent::computeQpJacobian()
{
  // The kernel variable does not appear in the residual (no _u[_qp]),
  // so the on-diagonal Jacobian is zero
  return 0.0;
}

Real
LorentzForceComponent::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Derivatives with respect to coupled variables
  // f = mu_0 * (J x H)
  // For component 0 (x): f_x = mu_0 * (Jy * Hz - Jz * Hy)
  //   dFx/dJy = mu_0 * Hz, dFx/dJz = -mu_0 * Hy
  //   dFx/dHy = -mu_0 * Jz, dFx/dHz = mu_0 * Jy
  // For component 1 (y): f_y = mu_0 * (Jz * Hx - Jx * Hz)
  //   dFy/dJz = mu_0 * Hx, dFy/dJx = -mu_0 * Hz
  //   dFy/dHx = mu_0 * Jz, dFy/dHz = -mu_0 * Jx
  // For component 2 (z): f_z = mu_0 * (Jx * Hy - Jy * Hx)
  //   dFz/dJx = mu_0 * Hy, dFz/dJy = -mu_0 * Hx
  //   dFz/dHx = -mu_0 * Jy, dFz/dHy = mu_0 * Jx

  if (_component == 0) // x-component: f_x = mu_0 * (Jy * Hz - Jz * Hy)
  {
    if (jvar == _Jy_var)
      return _mu_0 * _Hz[_qp] * _test[_i][_qp];
    else if (jvar == _Jz_var)
      return -_mu_0 * _Hy[_qp] * _test[_i][_qp];
    else if (jvar == _Hy_var)
      return -_mu_0 * _Jz[_qp] * _test[_i][_qp];
    else if (jvar == _Hz_var)
      return _mu_0 * _Jy[_qp] * _test[_i][_qp];
  }
  else if (_component == 1) // y-component: f_y = mu_0 * (Jz * Hx - Jx * Hz)
  {
    if (jvar == _Jz_var)
      return _mu_0 * _Hx[_qp] * _test[_i][_qp];
    else if (jvar == _Jx_var)
      return -_mu_0 * _Hz[_qp] * _test[_i][_qp];
    else if (jvar == _Hx_var)
      return _mu_0 * _Jz[_qp] * _test[_i][_qp];
    else if (jvar == _Hz_var)
      return -_mu_0 * _Jx[_qp] * _test[_i][_qp];
  }
  else if (_component == 2) // z-component: f_z = mu_0 * (Jx * Hy - Jy * Hx)
  {
    if (jvar == _Jx_var)
      return _mu_0 * _Hy[_qp] * _test[_i][_qp];
    else if (jvar == _Jy_var)
      return -_mu_0 * _Hx[_qp] * _test[_i][_qp];
    else if (jvar == _Hx_var)
      return -_mu_0 * _Jy[_qp] * _test[_i][_qp];
    else if (jvar == _Hy_var)
      return _mu_0 * _Jx[_qp] * _test[_i][_qp];
  }

  return 0.0;
}