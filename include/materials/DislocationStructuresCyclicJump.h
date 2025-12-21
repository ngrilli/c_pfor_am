// Marco Pelegatti
// Università di Udine
// Nicolò Grilli
// Università di Bristol
// 20 Dicembre 2025

#pragma once

#include "CrystalPlasticityCyclicDislocationStructures.h"

class DislocationStructuresCyclicJump;

/**
 * DislocationStructuresCyclicJump extrapolates variables
 * and allows for cyclic jumps to accelerate simulations
 */

class DislocationStructuresCyclicJump : public CrystalPlasticityCyclicDislocationStructures
{
public:
  static InputParameters validParams();

  DislocationStructuresCyclicJump(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
   */
  virtual void initQpStatefulProperties() override;

  virtual bool calculateSlipRate() override;

  const VariableValue & _cyclic_jump;  
};
