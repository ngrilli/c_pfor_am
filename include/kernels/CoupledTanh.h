// Nicolò Grilli
// Università di Bristol
// 10 Dicembre 2023

#pragma once

#include "Kernel.h"

// Hyperbolic tangent of a coupled variable

class CoupledTanh : public Kernel
{
public:
  static InputParameters validParams();

  CoupledTanh(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  
  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // The couple variable
  const VariableValue & _v;
  const bool _v_coupled;
  unsigned int _v_var;
  
  // Prefactor
  const Real _A;
  
  // Prefactor inside tanh
  const Real _theta;
  
  // Normalization constant for the coupled variable
  const Real _vn;
};
