// Nicolò Grilli
// University of Bristol
// 21 Luglio 2021

// DG upwinding for the advection of a coupled variable
// It is assumed that upwind scheme is always used
// Therefore, no need for a upwinding_type flag as in ConservativeAdvectionCoupled
// This kernel implements the term (on LHS): 
// d(rho_coupled v)/dx if dislo_character = edge
// d(rho_coupled v)/dy if dislo_character = screw
// Residual given by neighbouring coupled variable not considered

#pragma once

#include "DGKernel.h"

class DGAdvectionCoupledNN;

template <>
InputParameters validParams<DGAdvectionCoupledNN>();

class DGAdvectionCoupledNN : public DGKernel
{
public:
  static InputParameters validParams();

  DGAdvectionCoupledNN(const InputParameters & parameters);

protected:
  virtual void getDislocationVelocity();
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// advection velocity
  RealVectorValue _velocity;
  
  // Coupled dislocation density in the flux term
  const VariableValue & _rho_coupled;
  
  const bool _rho_coupled_coupled;
  unsigned int _rho_coupled_var;
  
  // Coupled dislocation density in the neighbouring element
  const VariableValue & _rho_neighbor;
  
  // Edge slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _edge_slip_direction;

  // Screw slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
  // Dislocation velocity value (signed) on all slip systems
  const MaterialProperty<std::vector<Real>> & _dislo_velocity;

  // Slip system index to determine slip direction
  const unsigned int _slip_sys_index;
  
  // Character of dislocations (edge or screw)
  const enum class DisloCharacter { edge, screw } _dislo_character;
  
  // Check if the variable is edge or screw, or total density
  // true = it is edge or screw
  // false = it is total density
  bool _is_edge_or_screw; 
  
};

