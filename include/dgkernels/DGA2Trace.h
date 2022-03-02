// Nicolo Grilli
// University of Bristol
// 20 Ottobre 2021

#pragma once

#include "DGKernel.h"

class DGA2Trace : public DGKernel
{
public:
  static InputParameters validParams();

  DGA2Trace(const InputParameters & parameters);

protected:
  virtual void getDislocationVelocity();
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// advection velocity direction
  /// not multiplied by velocity magnitude
  RealVectorValue _velocity;

  // Total dislocation density: rho_t 
  const VariableValue & _rho_tot;
  const bool _rho_tot_coupled;
  unsigned int _rho_tot_var;
  
  // Coupled dislocation density rho_t in the neighbouring element
  const VariableValue & _rho_neighbor;
  
  // Derivative of the velocity with respect to edge and screw slip directions
  const VariableValue & _dv_dx;
  const VariableValue & _dv_dy;
  
  // same in the neighbouring element 
  const VariableValue & _dv_dx_neighbor;
  const VariableValue & _dv_dy_neighbor;
  
  // Max absolute value of dv_dx and dv_dy
  const Real _dv_dx_max;
  const Real _dv_dy_max;  
  
  // Edge slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _edge_slip_direction;

  // Screw slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
  // Slip system index to determine slip direction
  const unsigned int _slip_sys_index;
  
  // Character of dislocations (edge or screw)
  const enum class DisloCharacter { edge, screw } _dislo_character;

};

