// Nicolò Grilli
// National University of Singapore
// 22 Novembre 2020

#pragma once

#include "Kernel.h"

class DisloMultiplication;

template <>
InputParameters validParams<DisloMultiplication>();

class DisloMultiplication : public Kernel
{
public:
  static InputParameters validParams();

  DisloMultiplication(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const VariableValue & _rho_mult_1;
  const VariableValue & _rho_mult_2;
  
  const bool _rho_mult_1_coupled;
  const bool _rho_mult_2_coupled;
  
  unsigned int _rho_mult_1_var;
  unsigned int _rho_mult_2_var;

};

