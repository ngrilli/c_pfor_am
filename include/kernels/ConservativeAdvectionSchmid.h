// Nicolo Grilli
// National University of Singapore
// 13 Ottobre 2020

#pragma once

#include "Kernel.h"

// Forward Declaration
class ConservativeAdvectionSchmid;

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 * Advection velocity \vec{v} is taken as material property.
 */
template <>
InputParameters validParams<ConservativeAdvectionSchmid>();

class ConservativeAdvectionSchmid : public Kernel
{
public:
  static InputParameters validParams();

  ConservativeAdvectionSchmid(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual void computeResidual() override;
  virtual void computeJacobian() override;

  /// advection velocity
  std::vector<Real> _velocity;

  /// enum to make the code clearer
  enum class JacRes
  {
    CALCULATE_RESIDUAL = 0,
    CALCULATE_JACOBIAN = 1
  };
  
  // Slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _slip_direction;

  /// Type of upwinding
  const enum class UpwindingType { none, full } _upwinding;
  
  // Slip system index to determine slip direction
  const unsigned int _slip_sys_index;

  /// Nodal value of u, used for full upwinding
  const VariableValue & _u_nodal;

  /// In the full-upwind scheme, whether a node is an upwind node
  std::vector<bool> _upwind_node;

  /// In the full-upwind scheme d(total_mass_out)/d(variable_at_node_i)
  std::vector<Real> _dtotal_mass_out;

  /// Returns - _grad_test * velocity
  Real negSpeedQp();

  /// Calculates the fully-upwind Residual and Jacobian (depending on res_or_jac)
  void fullUpwind(JacRes res_or_jac);
};
