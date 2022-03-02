// Nicolò Grilli
// University of Bristol
// 4 Settembre 2021

// DG upwinding for the advection of a coupled variable
// It is assumed that upwind scheme is always used
// Therefore, no need for a upwinding_type flag as in ConservativeAdvectionCoupled
// This kernel implements the term (on LHS): 
// d(rho_gnd v)/dx if dislo_character = edge
// d(rho_gnd v)/dy if dislo_character = screw
// Upwind condition is calculated both on edge/screw dislocations
// in this element and on the neighbouring element
// The forward and backward motion of
// positive and negative GND is taken into account.
// The velocity directions of GND with arbitrary character
// perpendicular to the line segment is taken into account.
// This kernel must be applied to rho_tot

#pragma once

#include "DGKernel.h"

class DGAdvectionRhoTot : public DGKernel
{
public:
  static InputParameters validParams();

  DGAdvectionRhoTot(const InputParameters & parameters);

protected:
  virtual void getDislocationVelocity();
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// advection velocity
  RealVectorValue _velocity;
  
  // Coupled edge dislocation density in the flux term
  const VariableValue & _rho_edge;
  
  const bool _rho_edge_coupled;
  unsigned int _rho_edge_var;
  
  // Coupled dislocation density in the neighbouring element
  const VariableValue & _rho_edge_neighbor;
  
  // Coupled screw dislocation density in the flux term
  // The following relationship applies for the variables in Hochrainer model
  // in the case of a purely signed dislocation (no SSD):
  // rho_t = sqrt(rho_e^2 + rho_s^2)
  const VariableValue & _rho_screw;
  
  const bool _rho_screw_coupled;
  unsigned int _rho_screw_var;
  
  // Coupled dislocation density of the other type in the neighbouring element
  const VariableValue & _rho_screw_neighbor;
  
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
  
  // Check that |rho_gnd| / rho_tot <= 1
  // In theory, total dislocation density cannot become higher than GND density
  // but it can happen because of numerical error
  bool _check_gnd_rho_ratio;
  
  // Tolerance on small values of rho_tot
  // This is needed for the velocity direction calculation
  const Real _rho_tot_tol;

};

