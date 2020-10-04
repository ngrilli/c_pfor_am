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
   * This function calculate stress residual.
   */
  virtual void calcResidual( RankTwoTensor &resid );

  const VariableValue & _temp;
  const Real _reference_temperature;
  const Real _thermal_expansion;
};


