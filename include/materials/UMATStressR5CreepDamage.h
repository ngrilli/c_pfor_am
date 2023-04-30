// Nicolò Grilli
// University of Bristol
// 30 Aprile 2023

#pragma once

#include "UMATStressDamage.h"

/**
 * Coupling material to use Abaqus UMAT models in MOOSE
 * with phase field damage
 * 
 * This is a variant of UMATStressDamage in which the fracture energy
 * is degraded based on the R5 creep damage criterion, which is based on
 * ductility exhaustion theory.
 * 
 * The damage model is reported in:
 * 
 * M. W. Spindler,
 * The prediction of creep damage in type 347 weld metal. 
 * Part I: the determination of material properties from creep and tensile tests
 * International Journal of Pressure Vessels and Piping
 * Volume 82, Issue 3, March 2005, Pages 175-184
 *
 */
class UMATStressR5CreepDamage : public UMATStressDamage
{
public:
  static InputParameters validParams();

  UMATStressR5CreepDamage(const InputParameters & parameters);

protected:

  /// Creep damage degradation function is added, which corresponds to a decrease of Gc
  /// in the damage model
  virtual void computeHistoryVariable(Real & F_pos, Real & F_neg);
  
  /// Minimum residual creep degradation
  const Real _residual_creep_degradation;
  
  /// creep degradation function
  const MaterialProperty<Real> & _creep_degradation_old;
  
  /// creep damage function
  /// calculated by the UMAT and passed to this MOOSE material property
  MaterialProperty<Real> & _f_ep_c;

};
