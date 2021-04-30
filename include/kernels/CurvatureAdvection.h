// Nicolo Grilli
// National University of Singapore
// 25 Gennaio 2021

#pragma once

#include "Kernel.h"

// Forward Declaration
class CurvatureAdvection;

/**
 * Dislocation curvature advection
 * first term on the right hand side of equation (3) in
 * Stefan Sandfeld and Michael Zaiser
 * Pattern formation in a minimal model of continuum dislocation plasticity
 * Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp)
 */
template <>
InputParameters validParams<CurvatureAdvection>();

class CurvatureAdvection : public Kernel
{
public:
  static InputParameters validParams();

  CurvatureAdvection(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// advection velocity
  std::vector<Real> _velocity;
  
  // GND dislocation density: rho_x or rho_y for edge or screw
  const VariableValue & _rho_gnd;
  const bool _rho_gnd_coupled;
  unsigned int _rho_gnd_var; 

  // Total dislocation density: rho_t 
  const VariableValue & _rho_tot;
  const bool _rho_tot_coupled;
  unsigned int _rho_tot_var;
  
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
  
  // Sign of dislocations
  const enum class DisloSign { positive, negative } _dislo_sign;
  
  // Character of dislocations (edge or screw)
  const enum class DisloCharacter { edge, screw } _dislo_character;

  /// Returns - _grad_test * velocity
  Real negSpeedQp();
};
