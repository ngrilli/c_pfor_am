// Nicolò Grilli
// University of Bristol
// 7 Giugno 2024

#include "DiscreteNucleationForceTemp.h"

registerMooseObject("c_pfor_amApp", DiscreteNucleationForceTemp);

InputParameters
DiscreteNucleationForceTemp::validParams()
{
  InputParameters params = DiscreteNucleationForce::validParams();
  params.addClassDescription(
      "Term for inserting grain nuclei or phases in non-conserved order parameter fields "
      "depending on the temperature. ");
  params.addCoupledVar("temperature", "Temperature. ");
  return params;
}

DiscreteNucleationForceTemp::DiscreteNucleationForceTemp(const InputParameters & params)
  : DiscreteNucleationForce(params),
  _temperature(coupledValue("temperature"))
{
}

Real
DiscreteNucleationForceTemp::computeQpResidual()
{
  // residual must be controlled by the temperature
  return -((*_nucleus)[_qp] * (_v1 - _v0) + _v0) * _test[_i][_qp];
}
