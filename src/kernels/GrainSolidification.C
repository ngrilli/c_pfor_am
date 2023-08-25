// Nicolò Grilli
// Università di Bristol
// 25 Agosto 2023

// Term representing the interaction between phases and grain orientations
// last term in equation (10) in:
// Min Yang, Lu Wang, Wentao Yan
// Phase-field modeling of grain evolutions in additive manufacturing from nucleation, growth, to coarsening
// https://www.nature.com/articles/s41524-021-00524-6

#include "GrainSolidification.h"

registerMooseObject("PhaseFieldApp", GrainSolidification);

InputParameters
GrainSolidification::validParams()
{
  InputParameters params = ACGrGrBase::validParams();
  params.addClassDescription("Solidification term in the grain free energy");
  params.addRequiredCoupledVar("zeta","Phase field representing liquid (0) or solid (1)");
  return params;
}

GrainSolidification::GrainSolidification(const InputParameters & parameters)
  : ACGrGrBase(parameters), 
  _gamma(getMaterialProperty<Real>("gamma_asymm")),
  _zeta(coupledValue("zeta")),
  _zeta_coupled(isCoupled("zeta")),
  _zeta_var(_zeta_coupled ? coupled("zeta") : 0)
{
}

Real
GrainSolidification::computeDFDOP(PFFunctionType type)
{
  // Calculate either the residual or Jacobian of the solidification term of the grain growth free energy
  switch (type)
  {
    case Residual:
    {
      return _mu[_qp] * (1.0 - _zeta[_qp]) * (1.0 - _zeta[_qp]) * 2.0 * _u[_qp];
    }

    case Jacobian:
    {
      return _mu[_qp] * (1.0 - _zeta[_qp]) * (1.0 - _zeta[_qp]) * 2.0 * _phi[_j][_qp];
    }

    default:
      mooseError("Invalid type passed in GrainSolidification");
  }
}

// Residual does not depend on the other phase field variables Etaj
Real
GrainSolidification::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real dDFDOP;

  if (_zeta_coupled && jvar == _zeta_var) {

    dDFDOP = (-4.0) * (1.0 - _zeta[_qp]) * _u[_qp] * _phi[_j][_qp];

    return _L[_qp] * _test[_i][_qp] * dDFDOP;
  }

  return 0.0;
}
