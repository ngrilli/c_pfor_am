// Nicolò Grilli
// Università di Bristol
// 14 Ottobre 2026

#pragma once

#include "Kernel.h"

class CurlCurlComponent : public Kernel
{
public:
  static InputParameters validParams();

  CurlCurlComponent(const InputParameters & parameters);

protected:
  virtual void computeCurlH(Real & CurlHx, Real & CurlHy, Real & CurlHz);
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  //virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Identity of the coupled variable
  /// whose derivatives have positive and negative sign in the curl component
  const unsigned int _H1_var;
  const unsigned int _H2_var;

  /// Gradient of the coupled variables
  const VariableGradient & _grad_H1;
  const VariableGradient & _grad_H2;

  /// Component of the curl vector to match
  const unsigned int _component;

  /// Name of the material property for current-dependent resistivity
  const MaterialProperty<Real> & _resistivity;
};
