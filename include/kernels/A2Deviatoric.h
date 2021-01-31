// Nicolo Grilli
// National University of Singapore
// 30 Gennaio 2021

#pragma once

#include "Kernel.h"

// Forward Declaration
class A2Deviatoric;

/**
 * Dislocation curvature diffusion
 * second term on the right hand side of equation (3) in
 * Stefan Sandfeld and Michael Zaiser (deviatoric part)
 * Pattern formation in a minimal model of continuum dislocation plasticity
 * Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp)
 */
template <>
InputParameters validParams<A2Deviatoric>();

class A2Deviatoric : public Kernel
{
public:
  static InputParameters validParams();

  A2Deviatoric(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// advection velocity
  std::vector<Real> _velocity;

  // Edge dislocation density: rho_x 
  const VariableValue & _rho_gnd_edge;
  const bool _rho_gnd_edge_coupled;
  unsigned int _rho_gnd_edge_var;
  
  // Screw dislocation density: rho_y 
  const VariableValue & _rho_gnd_screw;
  const bool _rho_gnd_screw_coupled;
  unsigned int _rho_gnd_screw_var;  
  
  // Derivative of the velocity with respect to edge and screw slip directions
  const VariableValue & _dv_dx;
  const VariableValue & _dv_dy;
  
  // Edge slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _edge_slip_direction;

  // Screw slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
  // Slip system index to determine slip direction
  const unsigned int _slip_sys_index;
  
  // Sign of dislocations
  const enum class DisloSign { positive, negative } _dislo_sign;
  
  // Character of dislocations (edge or screw)
  const enum class DisloCharacter { edge, screw } _dislo_character;

  /// Returns - _grad_test * velocity
  Real negSpeedQp();
};
