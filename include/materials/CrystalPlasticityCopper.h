// Nicolò Grilli
// Università di Bristol
// 2 Agosto 2025

#pragma once

#include "CrystalPlasticityKalidindiUpdate.h"

class CrystalPlasticityCopper;

/**
 * CrystalPlasticityCopper uses the multiplicative decomposition of the
 * deformation gradient and solves the PK2 stress residual equation at the
 * intermediate configuration to evolve the material state. The internal
 * variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Phenomenological model for copper with spatially random gss_initial.
 */

class CrystalPlasticityCopper : public CrystalPlasticityKalidindiUpdate
{
public:
  static InputParameters validParams();

  CrystalPlasticityCopper(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
   * Added random variations of the spatial distribution of the
   * initial lattice friction strength of the material
   */
  virtual void initQpStatefulProperties() override;

  /// Standard deviation of the spatial distribution of the
  /// initial lattice friction strength of the material
  const Real _gss_initial_std;
};
