// Nicolo Grilli
// University of Bristol
// Daijun Hu
// National University of Singapore
// 26 Luglio 2021

#pragma once

#include "DGKernel.h"

class DGCurvatureAdvection : public DGKernel
{
public:
  static InputParameters validParams();

  DGCurvatureAdvection(const InputParameters & parameters);

protected:
  virtual void getDislocationVelocity();
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// advection velocity
  RealVectorValue _velocity;
  
  // GND dislocation density: rho_x or rho_y for edge or screw
  const VariableValue & _rho_gnd;
  
  const bool _rho_gnd_coupled;
  unsigned int _rho_gnd_var; 

  // Total dislocation density: rho_t 
  const VariableValue & _rho_tot;
  
  const bool _rho_tot_coupled;
  unsigned int _rho_tot_var;
  
  // Coupled dislocation densities in the neighbouring element
  const VariableValue & _rho_gnd_neighbor;
  const VariableValue & _rho_tot_neighbor;

  // Tolerance on small values of rho_tot
  const Real _rho_tot_tol;
  
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
  
  // No need to check that variable _u on which the kernel is applied is
  // GND density or total density like in DGAdvectionCoupled
  
};
