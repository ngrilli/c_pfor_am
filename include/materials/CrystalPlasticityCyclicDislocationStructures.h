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
  
  /// Initialize constant interaction matrix between slip systems
  virtual void initializeInteractionMatrix();
  
  /// read Euler angles from file and calculate rotation matrix
  virtual void assignEulerAngles();

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
  
  // Backstress update based on (15) e (16)
  virtual void BackstressUpdate();

  /*
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual bool updateStateVariables() override;
  
  // Calculate wall volume fraction based on (44)
  virtual void calculateWallVolumeFraction();
  
  // Calculate substructure parameter eta based on (41), (42), (43)
  virtual void calculateSubstructureParameter();
  
  // Calculate substructure size based on (40)
  virtual void calculateSubstructureSize();
  
  // Calculate PSB fraction
  virtual void calculatePSBFraction();
  
  // Slip rate constants
  const Real _ao;
  const Real _xm;
  
  // Cap the absolute value of the slip increment in one time step to _slip_incr_tol
  const bool _cap_slip_increment;

  // Magnitude of the Burgers vector
  const Real _burgers_vector_mag;
  
  // Shear modulus in Taylor hardening law G
  const Real _shear_modulus;
  
  // Thermal slip resistance
  const Real _s_0;
  
  // Poisson's ratio for backstress calculation by Eshelby's inclusion
  const Real _nu;
  
  // Constant of similitude for dislocation substructure
  const Real _K_struct;
  
  // Resolved shear stress at initial yield
  const Real _tau_0;
  
  // Coefficient K in channel dislocations evolution, representing accumulation rate
  const Real _k_c;
  
  // Critical annihilation diameter for screw dislocations
  const Real _y_s;
  
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

  // Initial values of the channel dislocation density
  const Real _init_rho_c;
  
  // Initial values of the wall dislocation density
  const Real _init_rho_w;
  
  // Coefficient K in wall dislocations evolution, representing accumulation rate
  const Real _k_w;
  
  // Critical annihilation diameter for edge dislocations
  const Real _y_e;
  
  // Initial PSB fraction
  const Real _f_PSB_0;
  
  // PSB fraction at stabilization
  const Real _f_PSB_inf;
  
  // Increasing rate of PSB fraction
  const Real _k_PSB;
  
  // Critical accumulated plastic strain to develop PSBs
  const Real _epsilon_p_eff_cum_PSB;
  
  // Max/Min axis length ratio of the PSB
  const Real _eta_PSB;
  
  // PSB dislocation walls volume fraction
  const Real _f_w_PSB;
  
  // Initial characteristic PSB length
  const Real _init_d_struct_PSB;
  
  // Coefficient K in PSB dislocations evolution, representing accumulation rate
  const Real _k_c_PSB;
  
  // Critical annihilation diameter for dislocations in PSBs
  const Real _y_PSB;
  
  // Initial PSB dislocation density
  const Real _init_rho_PSB;
  
  // Dislocation densities: channel, wall, PSBs
  MaterialProperty<std::vector<Real>> & _rho_c;
  const MaterialProperty<std::vector<Real>> & _rho_c_old;
  MaterialProperty<std::vector<Real>> & _rho_w;
  const MaterialProperty<std::vector<Real>> & _rho_w_old;
  MaterialProperty<std::vector<Real>> & _rho_PSB;
  const MaterialProperty<std::vector<Real>> & _rho_PSB_old;
  
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
  const MaterialProperty<Real> & _d_struct_old;
  
  // Mean glide distance for dislocations in the channel phase
  MaterialProperty<Real> & _l_c;
  
  // Characteristic dislocation substructure length in the PSB
  MaterialProperty<Real> & _d_struct_PSB;
  const MaterialProperty<Real> & _d_struct_PSB_old;
  
  // Mean glide distance for dislocations in the PSB
  MaterialProperty<Real> & _l_PSB;
  
  // PSB fraction
  MaterialProperty<Real> & _f_PSB;
  const MaterialProperty<Real> & _f_PSB_old;
  
  // Backstress variables (channel, wall and PSB)
  MaterialProperty<std::vector<Real>> & _backstress_c;
  const MaterialProperty<std::vector<Real>> & _backstress_c_old;
  MaterialProperty<std::vector<Real>> & _backstress_w;
  const MaterialProperty<std::vector<Real>> & _backstress_w_old;
  MaterialProperty<std::vector<Real>> & _backstress_PSB;
  const MaterialProperty<std::vector<Real>> & _backstress_PSB_old;
  
  // Slip resistance in channel, wall and PSB
  MaterialProperty<std::vector<Real>> & _slip_resistance_c;
  const MaterialProperty<std::vector<Real>> & _slip_resistance_c_old;
  MaterialProperty<std::vector<Real>> & _slip_resistance_w;
  const MaterialProperty<std::vector<Real>> & _slip_resistance_w_old;
  MaterialProperty<std::vector<Real>> & _slip_resistance_PSB;
  const MaterialProperty<std::vector<Real>> & _slip_resistance_PSB_old;
  
  // Slip increment inside the channel, the wall and PSB
  MaterialProperty<std::vector<Real>> & _slip_increment_c;
  MaterialProperty<std::vector<Real>> & _slip_increment_w;
  MaterialProperty<std::vector<Real>> & _slip_increment_PSB;
  
  /// Increment of dislocation densities and backstress
  std::vector<Real> _rho_c_increment;
  std::vector<Real> _rho_w_increment;
  std::vector<Real> _rho_PSB_increment;
  std::vector<Real> _backstress_c_increment;
  std::vector<Real> _backstress_w_increment;
  std::vector<Real> _backstress_PSB_increment;

  /**
   * Stores the values of the dislocation densities and backstress 
   * from the previous substep
   * In classes which use dislocation densities, analogous dislocation density
   * substep vectors will be required.
   */
  std::vector<Real> _previous_substep_rho_c;
  std::vector<Real> _previous_substep_rho_w;
  std::vector<Real> _previous_substep_rho_PSB;
  std::vector<Real> _previous_substep_backstress_c;
  std::vector<Real> _previous_substep_backstress_w;
  std::vector<Real> _previous_substep_backstress_PSB;

  /**
   * Caches the value of the current dislocation densities immediately prior
   * to the update, and they are used to calculate the
   * the dislocation densities for the current substep (or step if
   * only one substep is taken) for the convergence check tolerance comparison.
   */
  std::vector<Real> _rho_c_before_update;
  std::vector<Real> _rho_w_before_update;
  std::vector<Real> _rho_PSB_before_update;
  std::vector<Real> _backstress_c_before_update;
  std::vector<Real> _backstress_w_before_update;
  std::vector<Real> _backstress_PSB_before_update;
  
  /// Interaction matrix between slip systems
  DenseMatrix<Real> _A_int;
  
  /// Intial macroscopic backstress: components and tensor
  std::vector<Real> _B_ii;
  RankTwoTensor _B_0;
  
  /// Element property read user object used to read initial dislocation structure size
  const PropertyReadFile * const _read_init_d;
  
  /// Element property read user object used to read in Euler angles
  const PropertyReadFile * const _read_prop_user_object;
  
  /// Values of the Euler Angles for rotation matrix calculation
  RealVectorValue _Euler_angles;
  
  /**
   * Crystal rotation in the original, or reference, configuration as defined by
   * Euler angle arguments
   */
  RankTwoTensor _crysrot;
  
};
