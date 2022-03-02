// Nicolò Grilli
// University of Bristol
// 7 Settembre 2021

// DG upwinding for the advection of a coupled variable
// It is assumed that upwind scheme is always used
// Therefore, no need for a upwinding_type flag as in ConservativeAdvectionCoupled
// This kernel implements the terms (on LHS): 
// d(rho_tot v)/dx if dislo_character = edge
// d(rho_tot v)/dy if dislo_character = screw
// Upwind condition is calculated both on edge/screw dislocations
// in this element and on the neighbouring element
// The forward and backward motion of
// positive and negative GND is taken into account.
// The velocity directions of GND with arbitrary character
// perpendicular to the line segment is taken into account.
// This kernel must be applied to rho_edge or rho_screw

#pragma once

#include "DGKernel.h"

class DGAdvectionRhoGND : public DGKernel
{
public:
  static InputParameters validParams();

  DGAdvectionRhoGND(const InputParameters & parameters);

protected:
  virtual void getDislocationVelocity();
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// advection velocity
  RealVectorValue _velocity;
  
  // Edge dislocation density
  const VariableValue & _rho_edge;
  
  const bool _rho_edge_coupled;
  unsigned int _rho_edge_var;
  
  // Edge dislocation density in the neighbouring element
  const VariableValue & _rho_edge_neighbor;
  
  // Screw dislocation density
  // The following relationship applies for the variables in Hochrainer model
  // in the case of a purely signed GND dislocation (no SSD):
  // rho_t = sqrt(rho_e^2 + rho_s^2)
  const VariableValue & _rho_screw;
  
  const bool _rho_screw_coupled;
  unsigned int _rho_screw_var;
  
  // Screw dislocation density in the neighbouring element
  const VariableValue & _rho_screw_neighbor;
  
  // Total dislocation density
  const VariableValue & _rho_tot;
  
  const bool _rho_tot_coupled;
  unsigned int _rho_tot_var;

  // Total dislocation density in the neighbouring element
  const VariableValue & _rho_tot_neighbor;
  
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

