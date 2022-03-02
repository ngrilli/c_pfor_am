// Nicolò Grilli
// University of Bristol
// 21 Luglio 2021

// DG upwinding for the advection of a coupled variable
// It is assumed that upwind scheme is always used
// Therefore, no need for a upwinding_type flag as in ConservativeAdvectionCoupled
// This kernel implements the term (on LHS): 
// d(rho_coupled v)/dx if dislo_character = edge
// d(rho_coupled v)/dy if dislo_character = screw
// Upwind condition is calculated both on edge/screw dislocations
// in this element and on the neighbouring element

#pragma once

#include "DGKernel.h"

class DGAdvectionCoupled : public DGKernel
{
public:
  static InputParameters validParams();

  DGAdvectionCoupled(const InputParameters & parameters);

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
  
  // Tolerance on small values of rho_tot
  // In this kernel it is used to avoid further outgoing dislocation flux
  // when the total dislocation density is very close to zero
  const Real _rho_tot_tol;
  
};

