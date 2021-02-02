// Nicolò Grilli
// National University of Singapore
// 2 Febbraio 2021

#pragma once

#include "Kernel.h"

/**
 * Residual of: coef * nabla^2 (coupled_variable)
 * After integration by part:
 * coef * (grad_test dot grad_some_var)
 */
class CoupledDiffusion : public Kernel
{
public:
  static InputParameters validParams();

  CoupledDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const VariableGradient & _velocity_vector;
  const bool _is_var_coupled;
  unsigned int _coupled_var;
  
  const Real & _coef;
};
