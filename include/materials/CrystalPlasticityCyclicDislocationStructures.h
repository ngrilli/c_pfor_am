// Marco Pelegatti
// Università di Udine
// Nicolò Grilli
// Università di Bristol
// 23 Marzo 2024

#pragma once

#include "CrystalPlasticityDislocationUpdateBase.h"
#include "ElementPropertyReadFile.h"

class CrystalPlasticityCyclicDislocationStructures;

/**
 * CrystalPlasticityCyclicDislocationStructures uses the multiplicative decomposition of the
 * deformation gradient and solves the PK2 stress residual equation at the
 * intermediate configuration to evolve the material state. The internal
 * variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Dislocation based model for crystal plasticity for cyclic load.
 * with backstress and dislocation structures based on:
 * Theodore Zirkle, Ting Zhu, David L. McDowell
 * Micromechanical crystal plasticity back stress evolution within FCC dislocation substructure
 * International Journal of Plasticity 146 (2021) 103082
 */

class CrystalPlasticityCyclicDislocationStructures : public CrystalPlasticityDislocationUpdateBase
{
public:
  static InputParameters validParams();

  CrystalPlasticityCyclicDislocationStructures(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
   */
  virtual void initQpStatefulProperties() override;

  /**
   * Sets the value of the current and previous substep iteration slip system
   * resistance to the old value at the start of the PK2 stress convergence
   * while loop.
   */
  virtual void setInitialConstitutiveVariableValues() override;

  /**
   * Sets the current slip system resistance value to the previous substep value.
   * In cases where only one substep is taken (or when the first) substep is taken,
   * this method just sets the current value to the old slip system resistance
   * value again.
   */
  virtual void setSubstepConstitutiveVariableValues() override;

  /**
   * Stores the current value of the slip system resistance into a separate
   * material property in case substepping is needed.
   */
  virtual void updateSubstepConstitutiveVariableValues() override;

  virtual bool calculateSlipRate() override;
  
  virtual void calculateSlipResistance();

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  // Cache the slip system value before the update for the diff in the convergence check
  virtual void cacheStateVariablesBeforeUpdate() override;

  // Evolution of the state variables 
  virtual void calculateStateVariableEvolutionRateComponent() override;
  
  // Armstrong-Frederick update of the backstress
  virtual void ArmstrongFrederickBackstressUpdate();

  /*
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual bool updateStateVariables() override;
  
  // Slip rate constants
  const Real _ao;
  const Real _xm;
  
  // Cap the absolute value of the slip increment in one time step to _slip_incr_tol
  const bool _cap_slip_increment;

  // Magnitude of the Burgers vector
  const Real _burgers_vector_mag;
  
  // Shear modulus in Taylor hardening law G
  const Real _shear_modulus;
  
  // Self hardening coefficient
  const Real _A_self;
  
  // Poisson's ratio for backstress calculation by Eshelby's inclusion
  const Real _nu;
  
  // Constant of similitude for dislocation substructure
  const Real _K_struct;
  
  // Resolved shear stress at initial yield
  const Real _tau_0;
  
  // Coefficient K in channel dislocations evolution, representing accumulation rate
  const Real _k_0;
  
  // Critical annihilation distance for screw dislocations
  const Real _y_c;
  
  // Initial dislocation walls volume fraction
  const Real _f_0;
  
  // Saturated dislocation walls volume fraction
  const Real _f_inf;
  
  // Rate constant for the evolution of the dislocation walls volume fraction
  const Real _k_f;
  
  // Initial Max/Min axis length ratio of the dislocation substructure
  const Real _eta_0;
  
  // Saturated Max/Min axis length ratio of the dislocation substructure
  const Real _eta_inf;
  
  // Normalization constant for the evolution of the dislocation substructure
  const Real _X_norm;
  
  // Initial characteristic dislocation substructure length
  const Real _init_d_struct;

  // Initial values of the dislocation density
  const Real _init_rho_c;
  
  // Dislocation densities
  MaterialProperty<std::vector<Real>> & _rho_c;
  const MaterialProperty<std::vector<Real>> & _rho_c_old;
  
  // Cumulative effective plastic strain
  const MaterialProperty<Real> & _epsilon_p_eff_cum;
  
  // Instantaneous plastic deformation tangent at the slip system level
  MaterialProperty<std::vector<Real>> & _dslip_dtau;
  
  // Dislocation walls volume fraction
  MaterialProperty<Real> & _f_w;
  
  // Max/min axis length ratio of the dislocation substructure
  MaterialProperty<Real> & _eta;
  
  // Characteristic dislocation substructure length
  MaterialProperty<Real> & _d_struct;
  
  // Mean glide distance for dislocations in the channel phase
  MaterialProperty<Real> & _l_c;
  
  // Backstress variables
  MaterialProperty<std::vector<Real>> & _backstress;
  const MaterialProperty<std::vector<Real>> & _backstress_old;
  
  /// Increment of dislocation densities and backstress
  std::vector<Real> _rho_c_increment;
  std::vector<Real> _backstress_increment;

  /**
   * Stores the values of the dislocation densities and backstress 
   * from the previous substep
   * In classes which use dislocation densities, analogous dislocation density
   * substep vectors will be required.
   */
  std::vector<Real> _previous_substep_rho_c;
  std::vector<Real> _previous_substep_backstress;

  /**
   * Caches the value of the current dislocation densities immediately prior
   * to the update, and they are used to calculate the
   * the dislocation densities for the current substep (or step if
   * only one substep is taken) for the convergence check tolerance comparison.
   */
  std::vector<Real> _rho_c_before_update;
  std::vector<Real> _backstress_before_update;
  
};
