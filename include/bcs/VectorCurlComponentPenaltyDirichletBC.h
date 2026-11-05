// Nicolò Grilli
// Università di Bristol
// 5 Novembre 2026

#pragma once

#include "IntegratedBC.h"

class Function;

/**
 * Enforces a penalty Dirichlet boundary condition for a specified component of the curl
 * of a vector field H = (Hx, Hy, Hz), where H is constructed from three coupled scalar
 * variables. The boundary condition penalizes the difference between the computed curl
 * component and a prescribed exact value.
 *
 * This is the component-wise similar of VectorCurlPenaltyDirichletBC for cases where
 * the vector field is represented by separate scalar variables rather than a single
 * vector variable.
 */
class VectorCurlComponentPenaltyDirichletBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  VectorCurlComponentPenaltyDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Helper to compute the specified component of curl(H) at the current quadrature point
  Real computeCurlComponent() const;

  /// The penalty coefficient
  const Real _penalty;

  /// The component of curl(H) to constrain (0, 1, 2)
  const unsigned int _component;

  /// Identity of the coupled variables representing Hx, Hy or Hz
  const unsigned int _H1_var;
  const unsigned int _H2_var;

  /// Gradients of the coupled variables
  const VariableGradient & _grad_H1;
  const VariableGradient & _grad_H2;

  /// The component functions for the exact curl value
  const Function & _exact_curl_function;
};