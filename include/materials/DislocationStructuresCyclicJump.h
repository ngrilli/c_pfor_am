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
  const VariableValue & _cyclic_jump_old;
  
  const Real _real_cycle_duration;
  
  MaterialProperty<std::vector<Real>> & _extrapolated_slip_increment_c;
  MaterialProperty<std::vector<Real>> & _extrapolated_slip_increment_w;
  MaterialProperty<std::vector<Real>> & _extrapolated_slip_increment_PSB;
  
  Real _time_begin_real_cycles;
};
