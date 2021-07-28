// Nicolò Grilli
// National University of Singapore
// 24 Gennaio 2021

// Last term on the right hand side of equation (1) in
// Stefan Sandfeld and Michael Zaiser
// Pattern formation in a minimal model of continuum
// dislocation plasticity
// Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp)

#pragma once

#include "Kernel.h"

class CurvatureMultiplication;

template <>
InputParameters validParams<CurvatureMultiplication>();

class CurvatureMultiplication : public Kernel
{
public:
  static InputParameters validParams();

  CurvatureMultiplication(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  // Opposite signed dislocation type leading to annihilation
  const VariableValue & _curvature;
  
  // Slip system index to determine slip system velocity
  const unsigned int _slip_sys_index;
  
  // Check that dislocation density is positive
  // if it went below zero, it should not be further decreased
  bool _check_rho_positive;  
  
  // Tolerance on small values of rho_tot
  // In this kernel it is used to avoid further decrease of the
  // total dislocation density below a threshold value
  const Real _rho_tot_tol;
  
  const bool _curvature_coupled;
  
  unsigned int _curvature_var;
  
  // Dislocation velocity value (signed) on all slip systems
  const MaterialProperty<std::vector<Real>> & _dislo_velocity;

};

