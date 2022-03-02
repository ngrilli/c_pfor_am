// Nicolo Grilli
// University of Bristol
// 6 Luglio 2021

#pragma once

#include "Kernel.h"

/**
 * Coupled advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */

class AdvectionCoupledVConst : public Kernel
{
public:
  static InputParameters validParams();

  AdvectionCoupledVConst(const InputParameters & parameters);

protected:
  /// Returns - _grad_test * velocity
  Real negSpeedQp() const;
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar) override;

  /// advection velocity
  const VectorVariableValue & _velocity;
  
  // Coupled dislocation density in the flux term
  const VariableValue & _rho_coupled;
  
  const bool _rho_coupled_coupled;
  unsigned int _rho_coupled_var; 

};
