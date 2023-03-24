// Nicolò Grilli
// University of Bristol
// 24 Marzo 2023

#pragma once

#include "ComputeCrystalPlasticityStressDamage.h"

//#include "CrystalPlasticityDislocationUpdateBase.h"
//#include "ComputeCrystalPlasticityEigenstrainBase.h"
//#include "ElementPropertyReadFile.h"

//#include "RankTwoTensor.h"
//#include "RankFourTensor.h"

/**
 * ComputeCrystalPlasticityCreepDamage (used together with CrystalPlasticityDislocationUpdateBase)
 * uses the multiplicative decomposition of the deformation gradient and solves the PK2 stress
 * residual equation at the intermediate configuration to evolve the material state. The internal
 * variables are updated using an iterative predictor-corrector algorithm. Backward Euler
 * integration rule is used for the rate equations.
 * Thermal eigenstrain is included
 *
 * The only difference between this class and ComputeMultipleCrystalPlasticityStress
 * is that the _models variable here is an array of CrystalPlasticityDislocationUpdateBase
 * instead of CrystalPlasticityStressUpdateBase.
 *
 * This material model is coupled with creep specific phase field damage
 *
 * Phase field damage formulation is based on:
 * Nicolo Grilli and Marisol Koslowski
 * The effect of crystal anisotropy and plastic response
 * on the dynamic fracture of energetic materials
 * Journal of Applied Physics 126, 155101 (2019).
 * Plastic work is added to the positive part of the free energy that causes damage
 * 
 * while creep damage formulation is based on:
 * 
 * Qikun Xie, Hongyu Qi, Shaolin Li, Xiaoguang Yang, Duoqi Shi, Fulin Li
 * Phase-field fracture modeling for creep crack
 * Theoretical and Applied Fracture Mechanics 124 (2023) 103798
 * 
 */
class ComputeCrystalPlasticityCreepDamage : public ComputeCrystalPlasticityStressDamage
{
public:
  static InputParameters validParams();

  ComputeCrystalPlasticityCreepDamage(const InputParameters & parameters);

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
   * This updates internal state variables once after the Newton-Raphson
   * crystal plasticity iterations
   */
  void postSolveQp(RankTwoTensor & stress_new, RankFourTensor & jacobian_mult);

  /// compute history variable and assign to _E
  /// which is used by the fracture model for damage growth
  /// Damage grows only because of the positive part of the elastic energy F_pos
  /// At this point the creep degradation function is introduced,
  /// this is equivalent to reducing the value of Gc									 
  virtual void computeHistoryVariable(Real & F_pos, Real & F_neg);
  
  /// Minimum residual creep degradation
  const Real _residual_creep_degradation;
  
  /// Creep degradation function that will be effectively
  /// a prefactor for Gc
  MaterialProperty<Real> & _creep_degradation;
  const MaterialProperty<Real> & _creep_degradation_old;

};
