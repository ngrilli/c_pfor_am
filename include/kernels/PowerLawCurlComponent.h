// Nicolò Grilli
// Università di Bristol
// 13 Settembre 2026

#pragma once

#include "Kernel.h"

/**
 * This kernel implements the curl(E) term in (1) of:
 * Bin Feng, Jun Ma, Ferdi van der Heiden, Zhixuan Zhang, Hanlin Zhu, Yameng Zhang, Jintao Hu and Nick Simpson
 * Three-dimensional electromagnetic–mechanical coupled modelling of multilayer HTS REBCO coated conductors considering strain effects
 * Superconductor Science and Technology, Volume 39, Number 3, 035008
 * https://iopscience.iop.org/article/10.1088/1361-6668/ae49bd 
 */

class PowerLawCurlComponent : public Kernel
{
public:
  static InputParameters validParams();

  PowerLawCurlComponent(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Identity of the coupled variables
  /// whose derivatives have positive and negative sign in the curl component
  const unsigned int _J1_var;
  const unsigned int _J2_var;

  /// Coupled variables and their gradients
  const VariableValue & _J1;
  const VariableValue & _J2;
  const VariableGradient & _grad_J1;
  const VariableGradient & _grad_J2;

  /// Component of the curl vector to match
  const unsigned int _component;

  /// Critical electric field
  const Real _E0;

  /// Critical current density
  const Real _J0;

  /// Power-law exponent
  const Real _n;
};
