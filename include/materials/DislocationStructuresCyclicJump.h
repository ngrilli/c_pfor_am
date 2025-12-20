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

  /**
   * Sets the current slip system resistance value to the previous substep value.
   * In cases where only one substep is taken (or when the first) substep is taken,
   * this method just sets the current value to the old slip system resistance
   * value again.
   */
  ///virtual void setSubstepConstitutiveVariableValues() override;

  /**
   * Stores the current value of the slip system resistance into a separate
   * material property in case substepping is needed.
   */
  ///virtual void updateSubstepConstitutiveVariableValues() override;

  ///virtual bool calculateSlipRate() override;
  
  ///virtual void calculateSlipResistance();

  ///virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  /// Cache the slip system value before the update for the diff in the convergence check
  ///virtual void cacheStateVariablesBeforeUpdate() override;

  /// Evolution of the state variables 
  ///virtual void calculateStateVariableEvolutionRateComponent() override;
  
  /// Backstress update based on (15) e (16)
  ///virtual void BackstressUpdate();

  /**
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  ///virtual bool updateStateVariables() override;
  
  /// Calculate wall volume fraction based on (44)
  ///virtual void calculateWallVolumeFraction();
  
  /// Calculate substructure parameter eta based on (41), (42), (43)
  ///virtual void calculateSubstructureParameter();
  
  /// Calculate substructure size based on (40)
  ///virtual void calculateSubstructureSize();
  
  /// Calculate PSB fraction
  ///virtual void calculatePSBFraction();
  
  /// Slip rate constants
  ///const Real _gamma_o;
  ///const Real _m_exp;
  
  /// Cap the absolute value of the slip increment in one time step to _slip_incr_tol
  ///const bool _cap_slip_increment;

  /// Magnitude of the Burgers vector
  ///const Real _burgers_vector_mag;
  
  /// Shear modulus in Taylor hardening law G
  ///const Real _shear_modulus;
  
  /// Initial slip resistance
  ///const Real _s_0;
  
  /// Initial slip resistance in the channel and wall
  ///const Real _s_0c;
  ///const Real _s_0w;
  
  /// Poisson's ratio for backstress calculation by Eshelby's inclusion
  ///const Real _nu;
  
  /// Coefficient of accumulation rate of dislocations in channel phase
  ///const Real _k_c;
  
  /// Critical annihilation distance for dislocations in channel phase
  ///const Real _y_c;
  
  /// Coefficient of accumulation rate of dislocations in wall phase
  ///const Real _k_w;
  
  /// Critical annihilation distance for dislocations in wall phase
  ///const Real _y_w;
  
  /// Initial wall volume fraction
  ///const Real _f_0;
  
  /// Saturated wall volume fraction
  ///const Real _f_inf;
  
  /// Rate constant for evolution of wall volume fraction
  ///const Real _k_f;
  
  /// Initial channel aspect ratio of dislocation structure
  ///const Real _eta_0;
  
  /// Saturated channel aspect ratio of dislocation structure
  ///const Real _eta_inf;
  
  /// Initial channel aspect ratio of cellular structure at slip system level
  ///const std::vector<Real> _eta_0_alpha;
  
  /// Normalization constant for evolution of dislocation structure
  ///const Real _X_norm;
  
  /// Constant of similitude for dislocation structure
  ///const Real _K_struct;
  
  /// Resolved shear stress at initial yield for similitude scaling law
  ///const Real _tau_0;
  
  /// Initial characteristic dislocation structure length
  ///const Real _init_d_struct;

  /// Initial channel dislocation density
  ///const Real _init_rho_c;
  
  /// Initial wall dislocation density
  ///const Real _init_rho_w;
  
  /// Initial PSB volume fraction
  ///const Real _f_PSB_0;
  
  /// Saturated PSB volume fraction
  ///const Real _f_PSB_inf;
  
  /// Rate constant for evolution of PSB volume fraction
  ///const Real _k_PSB;
  
  /// Critical cumulative effective plastic strain to develop PSB
  ///const Real _epsilon_p_eff_cum_PSB;
  
  /// Coefficient of accumulation rate of dislocations in PSB channel phase
  ///const Real _k_c_PSB;
  
  /// Critical annihilation distance for dislocations in PSB phase
  ///const Real _y_PSB;
  
  /// PSB wall volume fraction
  ///const Real _f_w_PSB;
  
  /// Channel aspect ratio of PSB
  ///const Real _eta_PSB;
  
  /// Initial characteristic PSB length
  ///const Real _init_d_struct_PSB;
  
  /// Initial PSB dislocation density
  ///const Real _init_rho_PSB;
  
  /// Interaction matrix coefficients between slip systems
  
  ///const Real _A_self;
  ///const Real _A_copl;
  ///const Real _A_CS;
  ///const Real _A_GJ;
  ///const Real _A_HL;
  ///const Real _A_LC;
  
  /// Dislocation densities: channel, wall, PSBs
  ///MaterialProperty<std::vector<Real>> & _rho_c;
  ///const MaterialProperty<std::vector<Real>> & _rho_c_old;
  ///MaterialProperty<std::vector<Real>> & _rho_w;
  ///const MaterialProperty<std::vector<Real>> & _rho_w_old;
  ///MaterialProperty<std::vector<Real>> & _rho_PSB;
  ///const MaterialProperty<std::vector<Real>> & _rho_PSB_old;
  
  /// Cumulative effective plastic strain
  ///const MaterialProperty<Real> & _epsilon_p_eff_cum;
  
  /// Instantaneous plastic deformation tangent at the slip system level
  ///MaterialProperty<std::vector<Real>> & _dslip_dtau;
  
  /// Wall volume fraction
  ///MaterialProperty<Real> & _f_w;
  
  /// Channel aspect ratio of dislocation structure: scalar or dependent on slip system
  ///MaterialProperty<Real> & _eta;
  ///MaterialProperty<std::vector<Real>> & _eta_vector;
  
  /// Characteristic dislocation structure length
  ///MaterialProperty<Real> & _d_struct;
  ///const MaterialProperty<Real> & _d_struct_old;
  
  /// Mean glide distance for dislocations in the channel phase: scalar or dependent on slip system
  ///MaterialProperty<Real> & _l_c;
  ///MaterialProperty<std::vector<Real>> & _l_c_vector;
  
  /// Characteristic PSB length
  ///MaterialProperty<Real> & _d_struct_PSB;
  ///const MaterialProperty<Real> & _d_struct_PSB_old;
  
  /// Mean glide distance for dislocations in the PSB
  ///MaterialProperty<Real> & _l_PSB;
  
  /// PSB volume fraction
  ///MaterialProperty<Real> & _f_PSB;
  ///const MaterialProperty<Real> & _f_PSB_old;
  
  /// Backstress variables in channel, wall and PSB
  ///MaterialProperty<std::vector<Real>> & _backstress_c;
  ///const MaterialProperty<std::vector<Real>> & _backstress_c_old;
  ///MaterialProperty<std::vector<Real>> & _backstress_w;
  ///const MaterialProperty<std::vector<Real>> & _backstress_w_old;
  ///MaterialProperty<std::vector<Real>> & _backstress_PSB;
  ///const MaterialProperty<std::vector<Real>> & _backstress_PSB_old;
  
  /// Slip resistance in channel, wall and PSB
  ///MaterialProperty<std::vector<Real>> & _slip_resistance_c;
  ///const MaterialProperty<std::vector<Real>> & _slip_resistance_c_old;
  ///MaterialProperty<std::vector<Real>> & _slip_resistance_w;
  ///const MaterialProperty<std::vector<Real>> & _slip_resistance_w_old;
  ///MaterialProperty<std::vector<Real>> & _slip_resistance_PSB;
  ///const MaterialProperty<std::vector<Real>> & _slip_resistance_PSB_old;
  
  /// Slip increment inside the channel, the wall and PSB
  ///MaterialProperty<std::vector<Real>> & _slip_increment_c;
  ///MaterialProperty<std::vector<Real>> & _slip_increment_w;
  ///MaterialProperty<std::vector<Real>> & _slip_increment_PSB;
  
  /// Increment of dislocation densities and backstress
  ///std::vector<Real> _rho_c_increment;
  ///std::vector<Real> _rho_w_increment;
  ///std::vector<Real> _rho_PSB_increment;
  ///std::vector<Real> _backstress_c_increment;
  ///std::vector<Real> _backstress_w_increment;
  ///std::vector<Real> _backstress_PSB_increment;

  /**
   * Stores the values of the dislocation densities and backstress 
   * from the previous substep
   * In classes which use dislocation densities, analogous dislocation density
   * substep vectors will be required.
   */
  ///std::vector<Real> _previous_substep_rho_c;
  ///std::vector<Real> _previous_substep_rho_w;
  ///std::vector<Real> _previous_substep_rho_PSB;
  ///std::vector<Real> _previous_substep_backstress_c;
  ///std::vector<Real> _previous_substep_backstress_w;
  ///std::vector<Real> _previous_substep_backstress_PSB;

  /**
   * Caches the value of the current dislocation densities immediately prior
   * to the update, and they are used to calculate the
   * the dislocation densities for the current substep (or step if
   * only one substep is taken) for the convergence check tolerance comparison.
   */
  ///std::vector<Real> _rho_c_before_update;
  ///std::vector<Real> _rho_w_before_update;
  ///std::vector<Real> _rho_PSB_before_update;
  ///std::vector<Real> _backstress_c_before_update;
  ///std::vector<Real> _backstress_w_before_update;
  ///std::vector<Real> _backstress_PSB_before_update;
  
  /// Interaction matrix between slip systems
  ///DenseMatrix<Real> _A_int;
  
  /// Intial macroscopic backstress: components and tensor
  ///std::vector<Real> _B_ii;
  ///RankTwoTensor _B_0;
  
  /// Element property read user object used to read initial dislocation structure size
  ///const PropertyReadFile * const _read_init_d;
  
  /// Element property read user object used to read in Euler angles
  ///const PropertyReadFile * const _read_prop_user_object;
  
  /// Values of the Euler Angles for rotation matrix calculation
  ///RealVectorValue _Euler_angles;
  
  /**
   * Crystal rotation in the original, or reference, configuration as defined by
   * Euler angle arguments
   */
  ///RankTwoTensor _crysrot;
  
};
