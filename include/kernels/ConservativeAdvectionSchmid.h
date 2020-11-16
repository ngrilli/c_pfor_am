// Nicolo Grilli
// National University of Singapore
// 16 Novembre 2020

#pragma once

#include "Kernel.h"

// Forward Declaration
class ConservativeAdvectionSchmid;

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 * Advection velocity direction \vec{v} and magnitude/sign are taken as material property.
 * Based on dislocation velocity model dependent on the resolved shear stress.
 * Signed edge and screw dislocations are considered
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
  
  // Edge slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _edge_slip_direction;

  // Screw slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
  // Dislocation velocity value (signed) on all slip systems
  const MaterialProperty<std::vector<Real>> & _dislo_velocity;

  /// Type of upwinding
  const enum class UpwindingType { none, full } _upwinding;
  
  // Slip system index to determine slip direction
  const unsigned int _slip_sys_index;
  
  // Sign of dislocations
  const enum class DisloSign { positive, negative } _dislo_sign;
  
  // Character of dislocations (edge or screw)
  const enum class DisloCharacter { edge, screw } _dislo_character;

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
