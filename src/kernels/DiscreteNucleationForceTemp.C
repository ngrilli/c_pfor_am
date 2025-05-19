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
      "depending on the zeta variable. ");
  params.addCoupledVar("zeta", "zeta variable which is 0 in the liquid phase and 1 in the solid phase. ");
  params.addParam<Real>("zeta_threshold", 0.9, "Threshold above which nucleation does not take place. ");
  return params;
}

DiscreteNucleationForceTemp::DiscreteNucleationForceTemp(const InputParameters & params)
  : DiscreteNucleationForce(params),
  _zeta(coupledValue("zeta")),
  _zeta_threshold(getParam<Real>("zeta_threshold"))
{
}

Real
DiscreteNucleationForceTemp::computeQpResidual()
{
  if (_zeta[_qp] < _zeta_threshold) {

    return -((*_nucleus)[_qp] * (_v1 - _v0) + _v0) * _test[_i][_qp];
    
  } else {
   
    return 0.0;
  }
}
