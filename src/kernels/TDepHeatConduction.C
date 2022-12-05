// Nicolò Grilli
// University of Bristol
// 5 Dicembre 2022

#include "TDepHeatConduction.h"
#include "MooseMesh.h"

registerMooseObjectAliased("HeatConductionApp", TDepHeatConduction);

InputParameters
TDepHeatConduction::validParams()
{
  InputParameters params = Diffusion::validParams();
  params.addClassDescription(
      "Computes residual/Jacobian contribution for $(k (T) \\nabla T, \\nabla \\psi)$ term. "
      "k(T) depends linearly on temperature. ");
  params.addParam<MaterialPropertyName>(
      "diffusion_coefficient",
      "thermal_conductivity",
      "Property name of the diffusivity (Default: thermal_conductivity)");
  params.addParam<Real>("dthermal_conductivity_dT",0.0,"Constant derivative of the thermal conductivity with respect to temperature. ");
  params.addParam<Real>("reference_temperature",293.0,"Reference temperature (K) at which thermal conductivity is _diffusion_coefficient[_qp]. ");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

TDepHeatConduction::TDepHeatConduction(const InputParameters & parameters)
  : Diffusion(parameters),
    _diffusion_coefficient(getMaterialProperty<Real>("diffusion_coefficient")),
    _dthermal_conductivity_dT(getParam<Real>("dthermal_conductivity_dT")),
    _reference_temperature(getParam<Real>("reference_temperature"))
{
}

Real
TDepHeatConduction::computeQpResidual()
{
  return (_diffusion_coefficient[_qp] + _dthermal_conductivity_dT * (_u[_qp] - _reference_temperature)) * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
TDepHeatConduction::computeQpJacobian()
{
  Real TDepK = _diffusion_coefficient[_qp] + _dthermal_conductivity_dT * (_u[_qp] - _reference_temperature);	
  
  Real jac = _dthermal_conductivity_dT * _phi[_j][_qp] * _grad_u[_qp] * _grad_test[_i][_qp]
           + TDepK * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  return jac;
}
