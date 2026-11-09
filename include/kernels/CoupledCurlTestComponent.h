// Nicolò Grilli
// Università di Bristol
// 9 Novembre 2026

#pragma once

#include "Kernel.h"

class CoupledCurlTestComponent : public Kernel
{
public:
  static InputParameters validParams();

  CoupledCurlTestComponent(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Identity of the coupled variables
  /// which are the specified components of the vector whose curl is calculated
  const unsigned int _Hx_var;
  const unsigned int _Hy_var;
  const unsigned int _Hz_var;

  /// Values of the coupled variables
  const VariableValue & _Hx;
  const VariableValue & _Hy;
  const VariableValue & _Hz;

  /// Component of H * curl(test) that is calculated by the kernel (0, 1, 2)
  const unsigned int _component;
};