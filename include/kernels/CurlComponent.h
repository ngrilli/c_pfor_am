// Nicolò Grilli
// Università di Bristol
// 12 Settembre 2026

#pragma once

#include "Kernel.h"

class CurlComponent : public Kernel
{
public:
  static InputParameters validParams();

  CurlComponent(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Identity of the coupled variable
  /// whose derivatives have positive and negative sign in the curl component
  const unsigned int _H1_var;
  const unsigned int _H2_var;

  /// Gradient of the coupled variables
  const VariableGradient & _grad_H1;
  const VariableGradient & _grad_H2;

  /// Component of the curl vector to match
  const unsigned int _component;
};
