// Nicolò Grilli
// Università di Bristol
// 4 Novembre 2026

#pragma once

#include "Kernel.h"

/**
 * Calculates one component of the magnetic Lorentz force density:
 * f = mu_0 * (J x H)
 * where J is the current density vector and H is the magnetic field strength vector.
 * The force is per unit volume.
 */
class LorentzForceComponent : public Kernel
{
public:
  static InputParameters validParams();

  LorentzForceComponent(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // Coupled variables: magnetic field components
  const unsigned int _Hx_var;
  const unsigned int _Hy_var;
  const unsigned int _Hz_var;

  // Coupled variables: current density components
  const unsigned int _Jx_var;
  const unsigned int _Jy_var;
  const unsigned int _Jz_var;

  // Values of the coupled variables at the quadrature point
  const VariableValue & _Hx;
  const VariableValue & _Hy;
  const VariableValue & _Hz;
  const VariableValue & _Jx;
  const VariableValue & _Jy;
  const VariableValue & _Jz;

  // Component of the Lorentz force vector to compute (0=x, 1=y, 2=z)
  const unsigned int _component;

  // Magnetic permeability (mu_0)
  const Real _mu_0;
};