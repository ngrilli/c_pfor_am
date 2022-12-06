// Nicolò Grilli
// University of Bristol
// 6 Dicembre 2022

#pragma once

#include "HyperElasticPhaseFieldIsoDamage.h"

/**
 * This class solves visco plastic model based on isotropically damaged stress
 * The damage parameter is obtained from phase field fracture kernel
 * Computes undamaged elastic strain energy and associated tensors used in phase field fracture
 * kernel
 * Thermal stress is included.
 */
class HyperElasticThermalDamage : public HyperElasticPhaseFieldIsoDamage
{
public:
  static InputParameters validParams();

  HyperElasticThermalDamage(const InputParameters & parameters);

protected:
  /**
   * This function computes PK2 stress modified to account for damage
   * Computes numerical stiffness if flag is true
   * Computes undamaged elastic strain energy and associated tensors
   * Includes thermal stress.
   */
  virtual void computeDamageStress();

  /// Reference temperature at which thermal expansion is zero.
  Real _reference_temperature;
  
  /// Thermal expansion coefficient.
  Real _thermal_expansion_coeff;

  /// Coupled temperature
  const VariableValue & _temperature;

};
