// Nicolo Grilli
// Daijun Hu 
// National University of Singapore
// 27 Ottobre 2020

#pragma once

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticity uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */
class FiniteStrainCrystalPlasticityThermal : public FiniteStrainCrystalPlasticity
{
public:
  static InputParameters validParams();

  FiniteStrainCrystalPlasticityThermal(const InputParameters & parameters);

protected:
  /**
   * This function calculates stress residual.
   */
  virtual void calcResidual( RankTwoTensor &resid );
  
  /**
  * This function updates the slip increments.
  * And derivative of slip w.r.t. resolved shear stress.
  */
  virtual void getSlipIncrements();
  
  /**
  * This function
  * stores the dislocation velocity
  * to couple with dislocation transport
  */
  virtual void OutputSlipDirection();

  const VariableValue & _temp;
  const Real _thermal_expansion;
  const Real _reference_temperature;

  const Real _dCRSS_dT_A;
  const Real _dCRSS_dT_B;
  const Real _dCRSS_dT_C;
  const Real _dCTE_dT;
  
  // critical resolved shear stress
  // exponentially decreased with temperature
  std::vector<Real> _gssT;
  
  // Green-Lagrange strain tensor expressed in the
  // lattice coordinate system
  MaterialProperty<RankTwoTensor> & _lattice_strain; 
  
  // Rotated slip direction to couple with dislocation transport
  // to indicate dislocation velocity direction for all slip systems
  MaterialProperty<std::vector<Real>> & _slip_direction;
  
  // Slip increment for output
  MaterialProperty<std::vector<Real>> & _slip_incr_out;

};


