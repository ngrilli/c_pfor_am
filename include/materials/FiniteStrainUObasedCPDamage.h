// Nicolò Grilli
// University of Bristol
// 26 Settembre 2021

#pragma once

#include "FiniteStrainUObasedCP.h"

// Are these needed?
#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "CrystalPlasticityStateVariable.h"
#include "CrystalPlasticityStateVarRateComponent.h"

/**
 * FiniteStrainUObasedCPDamage is the coupling between
 * crystal plasticity and damage
 */
class FiniteStrainUObasedCPDamage : public FiniteStrainUObasedCP
{
public:
  static InputParameters validParams();

  FiniteStrainUObasedCPDamage(const InputParameters & parameters);

protected:
  /**
   * calculate stress residual.
   * Damage is added to the stress calculation
   */
  virtual void calcResidual();

  /// Variable defining the phase field damage parameter
  const VariableValue & _c;

};
