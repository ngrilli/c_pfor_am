// Nicolò Grilli
// University of Bristol
// 7 Giugno 2024

#pragma once

#include "DiscreteNucleationForce.h"

/**
 * Free energy penalty contribution to force the nucleation of subresolution particles
 */
class DiscreteNucleationForceTemp : public DiscreteNucleationForce
{
public:
  static InputParameters validParams();

  DiscreteNucleationForceTemp(const InputParameters & params);

  Real computeQpResidual() override;

protected:

  const VariableValue & _temperature;
  
};
