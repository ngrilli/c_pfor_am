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
  
  /**
   * Sets the value of the current and previous substep iteration slip system
   * resistance to the old value at the start of the PK2 stress convergence
   * while loop.
   */
  virtual void setInitialConstitutiveVariableValues() override;
  
  /**
   * Sets the current state variables value to the previous substep value.
   * In cases where only one substep is taken (or when the first) substep is taken,
   * this method just sets the current value to the old state variables
   * value again.
   */
  virtual void setSubstepConstitutiveVariableValues() override;

  /// Standard deviation of the spatial distribution of the
  /// initial lattice friction strength of the material
  const Real _gss_initial_std;
  
  /// Backstress parameters
  const Real _h;
  const Real _h_D;
  
  /// Backstress variables
  MaterialProperty<std::vector<Real>> & _backstress;
  const MaterialProperty<std::vector<Real>> & _backstress_old;
  
  /// Increment of state variables
  std::vector<Real> _backstress_increment;

  /**
   * Stores the values of the state variables
   * from the previous substep
   */
  std::vector<Real> _previous_substep_backstress;

  /**
   * Caches the value of the current state variables immediately prior
   * to the update, and they are used to calculate the
   * the state variables for the current substep (or step if
   * only one substep is taken) for the convergence check tolerance comparison
   */
  std::vector<Real> _backstress_before_update;
  
};
