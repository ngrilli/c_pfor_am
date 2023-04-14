// Nicolò Grilli
// Edward Horton
// University of Bristol
// 12 Aprile 2023

#pragma once

#include "ComputeCrystalPlasticityStressDamage.h"

/**
 * This is a variant of ComputeCrystalPlasticityStressDamage in which the fracture energy
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
class ComputeCPStressR5CreepDamage : public ComputeCrystalPlasticityStressDamage
{
public:
  static InputParameters validParams();

  ComputeCPStressR5CreepDamage(const InputParameters & parameters);

protected:
  
  /**
   * initializes the stateful properties such as PK2 stress, resolved shear
   * stress, plastic deformation gradient, slip system resistances, 
   * creep degradation function, etc.
   * This class is often overwritten by inherting classes.
   */
  virtual void initQpStatefulProperties() override;
  
  /**
   * Save the final stress and internal variable values after the iterative solve.
   */
  void postSolveQp(RankTwoTensor & stress_new, RankFourTensor & jacobian_mult);
					   
  /// compute history variable and assign to _E
  /// which is used by the fracture model for damage growth
  /// Damage grows only because of the positive part of the elastic energy F_pos									 
  virtual void computeHistoryVariable(Real & F_pos, Real & F_neg);

  /// Minimum residual creep degradation
  const Real _residual_creep_degradation;
  
  /// Creep degradation function that will be effectively
  /// a prefactor for Gc
  MaterialProperty<Real> & _creep_degradation;
  const MaterialProperty<Real> & _creep_degradation_old;

};
