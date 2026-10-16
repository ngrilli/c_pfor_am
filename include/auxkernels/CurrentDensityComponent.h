// Nicolò Grilli
// Università di Bristol
// 16 Ottobre 2026

#pragma once

#include "AuxKernel.h"

/**
 * Calculates one component of the current density:
 * J = curl(H)
 */
class CurrentDensityComponent : public AuxKernel
{
public:
  static InputParameters validParams();

  CurrentDensityComponent(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  // Gradients of the magnetic field components
  const VariableGradient & _grad_Hx;
  const VariableGradient & _grad_Hy;
  const VariableGradient & _grad_Hz;

  // Component of the curl vector to match
  const unsigned int _component;
};
