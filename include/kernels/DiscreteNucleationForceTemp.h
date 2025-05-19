// Nicolò Grilli
// University of Bristol
// 7 Giugno 2024

#pragma once

#include "DiscreteNucleationForce.h"

/**
 * Free energy penalty contribution to force the nucleation of subresolution particles.
 * Depedence on zeta variable to avoid grain nucleation in the solid.
 */
class DiscreteNucleationForceTemp : public DiscreteNucleationForce
{
public:
  static InputParameters validParams();

  DiscreteNucleationForceTemp(const InputParameters & params);

  Real computeQpResidual() override;

protected:

  // zeta variable which is 0 in the liquid phase and 1 in the solid phase
  const VariableValue & _zeta;
  
  // Threshold above which nucleation does not take place
  const Real _zeta_threshold;
};
