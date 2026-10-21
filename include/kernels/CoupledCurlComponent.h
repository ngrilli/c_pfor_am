// Nicolò Grilli
// Università di Bristol
// 21 Ottobre 2026

#pragma once

#include "Kernel.h"

class CoupledCurlComponent : public Kernel
{
public:
  static InputParameters validParams();

  CoupledCurlComponent(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Identity of the coupled variables
  /// which are the specified components of the vector whose curl is calculated
  const unsigned int _Hx_var;
  const unsigned int _Hy_var;
  const unsigned int _Hz_var;

  /// Gradient of the coupled variables
  const VariableGradient & _grad_Hx;
  const VariableGradient & _grad_Hy;
  const VariableGradient & _grad_Hz;

  /// Component of the curl vector to match
  const unsigned int _component;
};
