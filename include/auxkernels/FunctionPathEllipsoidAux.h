// Nicolò Grilli
// University of Bristol
// 7 Agosto 2022

#pragma once

#include "AuxKernel.h"

/**
 * This AuxKernel increases an AuxVariable from 0 to 1 if the qp is inside
 * an ellipsoid that is moving according to a path provided through a function.
 * It can be applied to a level set variable
 * to simulate the material deposition during wire arc additive manufacturing
 * together with ActDeactElementsCoupled.
 */
class FunctionPathEllipsoidAux : public AuxKernel
{
public:
  static InputParameters validParams();

  FunctionPathEllipsoidAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  Real getIntegralValue();

  std::vector<const VariableValue *> _coupled_vars;
  Real _coef;
  unsigned int _order;
  std::vector<Real> _integration_coef;

  /// The old variable value (zero if order == 3)
  const VariableValue & _u_old;
  /// The older variable value (zero if order != 3)
  const VariableValue & _u_older;
};
