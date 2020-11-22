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
  
  // Slip system index to determine slip system velocity
  const unsigned int _slip_sys_index;
  
  // Characteristic length of dislocation loops
  const Real _Lc;
  
  const bool _rho_mult_1_coupled;
  const bool _rho_mult_2_coupled;
  
  unsigned int _rho_mult_1_var;
  unsigned int _rho_mult_2_var;
  
  // Dislocation velocity value (signed) on all slip systems
  const MaterialProperty<std::vector<Real>> & _dislo_velocity;

};

