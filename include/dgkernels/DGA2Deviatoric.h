// Nicolo Grilli
// University of Bristol
// 22 Ottobre 2021

#pragma once

#include "DGKernel.h"

class DGA2Deviatoric : public DGKernel
{
public:
  static InputParameters validParams();

  DGA2Deviatoric(const InputParameters & parameters);

protected:
  virtual void getDislocationVelocity();
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// advection velocity direction
  /// not multiplied by velocity magnitude
  RealVectorValue _velocity;

  // Edge dislocation density: rho_x
  const VariableValue & _rho_gnd_edge;
  const bool _rho_gnd_edge_coupled;
  unsigned int _rho_gnd_edge_var;
  
  // Edge dislocation density rho_x in the neighbouring element
  const VariableValue & _rho_gnd_edge_neighbor;
  
  // Screw dislocation density: rho_y
  const VariableValue & _rho_gnd_screw;
  const bool _rho_gnd_screw_coupled;
  unsigned int _rho_gnd_screw_var;  
  
  // Screw dislocation density rho_y in the neighbouring element
  const VariableValue & _rho_gnd_screw_neighbor;
  
  // Derivative of the velocity with respect to edge and screw slip directions
  const VariableValue & _dv_dx;
  const VariableValue & _dv_dy;
  
  // same in the neighbouring element 
  const VariableValue & _dv_dx_neighbor;
  const VariableValue & _dv_dy_neighbor;
  
  // Max absolute value of dv_dx and dv_dy
  const Real _dv_dx_max;
  const Real _dv_dy_max;  
  
  // Tolerance on small values of ksabs
  const Real _ksabs_tol;
  
  // Edge slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _edge_slip_direction;

  // Screw slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
  // Slip system index to determine slip direction
  const unsigned int _slip_sys_index;
  
  // Character of dislocations (edge or screw)
  const enum class DisloCharacter { edge, screw } _dislo_character;

};

