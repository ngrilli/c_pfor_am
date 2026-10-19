// Nicolò Grilli
// Università di Bristol
// 16 Ottobre 2026

#include "CurrentDensityComponent.h"

registerMooseObject("c_pfor_amApp", CurrentDensityComponent);

InputParameters
CurrentDensityComponent::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates one component of the current density: J = curl(H)");
  params.addRequiredCoupledVar("Hx","x-component of the magnetic field");
  params.addRequiredCoupledVar("Hy","y-component of the magnetic field");
  params.addRequiredCoupledVar("Hz","z-component of the magnetic field");
  params.addRequiredParam<unsigned int>("component","Current density component: 0 = Jx, 1 = Jy, 2 = Jz");
  return params;
}

CurrentDensityComponent::CurrentDensityComponent(const InputParameters & parameters)
  : AuxKernel(parameters),
    _grad_Hx(coupledGradient("Hx")),
    _grad_Hy(coupledGradient("Hy")),
    _grad_Hz(coupledGradient("Hz")),
    _component(getParam<unsigned int>("component"))
{
  if (_component > 2)
    paramError("component", "Must be 0 (Jx), 1 (Jy), or 2 (Jz).");
}

Real
CurrentDensityComponent::computeValue()
{
  if (_component == 0)
    return _grad_Hz[_qp](1) - _grad_Hy[_qp](2);
  else if (_component == 1)
    return _grad_Hx[_qp](2) - _grad_Hz[_qp](0);
  else if (_component == 2)
    return _grad_Hy[_qp](0) - _grad_Hx[_qp](1);
  else
    return 0.0; // This case should never happen due to the check in the constructor
}