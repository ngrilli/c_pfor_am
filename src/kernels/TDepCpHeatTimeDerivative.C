// Nicolò Grilli
// University of Bristol
// 4 Dicembre 2022

#include "TDepCpHeatTimeDerivative.h"

registerMooseObject("HeatConductionApp", TDepCpHeatTimeDerivative);

InputParameters
TDepCpHeatTimeDerivative::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription(
      "Time derivative term $\\rho c_p (T) \\frac{\\partial T}{\\partial t}$ of "
      "the heat equation for specific heat $c_p$ that depends linearly on temperature "
      "and the density $\\rho$.");

  // Density may be changing with deformation, so we must integrate
  // over current volume by setting the use_displaced_mesh flag.
  params.set<bool>("use_displaced_mesh") = true;

  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property "
                                        "at the reference temperature. ");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  params.addParam<Real>("dspecific_heat_dT",0.0,"Constant derivative of the specific heat with respect to temperature. ");
  params.addParam<Real>("reference_temperature",293.0,"Reference temperature (K) at which specific heat is _specific_heat[_qp]. ");
  return params;
}

TDepCpHeatTimeDerivative::TDepCpHeatTimeDerivative(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density_name")),
    _dspecific_heat_dT(getParam<Real>("dspecific_heat_dT")),
    _reference_temperature(getParam<Real>("reference_temperature"))
{
}

Real
TDepCpHeatTimeDerivative::computeQpResidual()
{
  return _density[_qp] * (_specific_heat[_qp] + _dspecific_heat_dT * (_u[_qp] - _reference_temperature)) * _test[_i][_qp] * _u_dot[_qp];
}

Real
TDepCpHeatTimeDerivative::computeQpJacobian()
{
  return _specific_heat[_qp] * _density[_qp] * TimeDerivative::computeQpJacobian();
}
