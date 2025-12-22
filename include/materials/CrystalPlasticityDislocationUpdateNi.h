// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// University of Bristol
// 04 November 2025

#pragma once

#include "CrystalPlasticityDislocationUpdateBase.h"
#include "ElementPropertyReadFile.h"

class CrystalPlasticityDislocationUpdateNi;

/**
 * CrystalPlasticityDislocationUpdate uses the multiplicative decomposition of the
 * deformation gradient and solves the PK2 stress residual equation at the
 * intermediate configuration to evolve the material state. The internal
 * variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Dislocation based model for crystal plasticity.
 * Slip, creep and backstress are included 
 */

class CrystalPlasticityDislocationUpdateNi : public CrystalPlasticityDislocationUpdateBase
{
public:
  static InputParameters validParams();

  CrystalPlasticityDislocationUpdateNi(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
   */
  virtual void initQpStatefulProperties() override;
  
  /**
   * A helper method to rotate the a direction and plane normal system set into
   * the local crystal lattice orientation as defined by the crystal rotation
   * tensor from the Elasticity tensor class.
   * Edge and screw slip directions are also assigned
   */
  virtual void calculateSchmidTensor(const unsigned int & number_dislocation_systems,
                             const std::vector<RealVectorValue> & plane_normal_vector,
                             const std::vector<RealVectorValue> & direction_vector,
                             std::vector<RankTwoTensor> & schmid_tensor,
                             const RankTwoTensor & crysrot);


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

  virtual void
  calculateEquivalentSlipIncrement(RankTwoTensor & /*equivalent_slip_increment*/) override;

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  // Cache the slip system value before the update for the diff in the convergence check
  virtual void cacheStateVariablesBeforeUpdate() override;

  /**
   * Following the Constitutive model for slip system resistance as given in
   * Dislocation based model:
   * E. Demir, I Gutierrez-Urrutia
   * Investigation of strain hardening near grain
   * boundaries of an aluminium oligocrystal:
   * Experiments and crystal based finite element method
   * International Journal of Plasticity 136 (2021) 102898
   */
  virtual void calculateStateVariableEvolutionRateComponent() override;
  
  // Armstrong-Frederick update of the backstress
  virtual void ArmstrongFrederickBackstressUpdate();

  /*
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual bool updateStateVariables() override;

  /*
   * Determines if the state variables, e.g. defect densities, have converged
   * by comparing the change in the values over the iteration period.
   */
  virtual bool areConstitutiveStateVariablesConverged() override;

  // Variables used in
  // Eralp Demir, Ivan Gutierrez-Urrutia
  // Investigation of strain hardening near grain boundaries of an aluminum oligocrystal: 
  // Experiments and crystal based finite element method
  // International Journal of Plasticity
  // Volume 136, January 2021, 102898
  
  // Slip rate constants
  const Real _ao;
  const Real _xm;
  
  // Optional xm material property for exponent for slip rate
  const bool _include_xm_matprop;
  const MaterialProperty<Real> * _xm_matprop;
  
  // Optional function for slip prefactor. If provided, the slip prefactor can be set as a function of time
  // This is useful for an initial plastic deformation followed by creep load
  const Function * const _ao_function;
  
  // Use Kocks 1976 temperature dependence for xm
  // this is equation (19) in
  // E.D. Cyr et al.
  // A three dimensional (3D) thermo-elasto-viscoplastic constitutive model for FCC polycrystals
  // International Journal of Plasticity 70 (2015) 166-190
  // https://www.sciencedirect.com/science/article/pii/S0749641915000571
  bool _use_kocks_T_dependence_for_xm;
  
  // Activate creep strain rate
  bool _creep_activated;
  
  // Creep rate constants
  const Real _creep_ao;
  const Real _creep_xm;
  
  // Optional function for creep prefactor. If provided, the creep prefactor can be set as a function of time
  // This is useful for an initial plastic deformation followed by creep load
  const Function * const _creep_ao_function;
  
  // Optional function for creep resistance. If provided, the creep resistance can be set as a function of time
  // This is useful for differentiating resistance for slip and creep
  const Function * const _creep_resistance_function;
  
  // Tertiary creep constants
  const Real _m_exponent;
  const Real _creep_t0;
  const Real _creep_t_denominator;
  
  // Cap the absolute value of the slip increment in one time step to _slip_incr_tol
  const bool _cap_slip_increment;

  // Magnitude of the Burgers vector
  const Real _burgers_vector_mag;
  
  // Shear modulus in Taylor hardening law G
  const Real _shear_modulus;
  const Real _dshear_dT;  // NEW
  
  // Prefactor of Taylor hardening law, alpha
  const Real _alpha_0;
  
  // Latent hardening coefficient
//  const Real _r;
  const Real _r_self; // NEW
  const Real _r_oct_oct; // NEW
  const Real _r_cub_cub; // NEW
  const Real _r_oct_cub; // NEW

// Peierls stress (split by slip type)
  const Real _tau_c_0_oct; // NEW: Peierls stress for {111} slip
  const Real _tau_c_0_cub; // NEW: Peierls stress for {100} slip
  
  // Optional function for Peierls stress. If provided, the Peierls stress can be set as a function of time
  // This is useful for time dependent solid solution strenthening and precipitation hardening, e.g. during thermal ageing
  const Function * const _tau_c_0_function;

  // Solid solution strengthening contribution (new for IN718)
  const Real _tau_ss; 
  
  // Hall-Petch strengthening parameters (new for IN718)
  const Real _k_hp; // <-- Hall-Petch
  // const Real _d_grain; 
  
  // Coefficient K in SSD evolution, representing accumulation rate
  const Real _k_0;
  
  // Critical annihilation diameter
  const Real _y_c;
  const Real _Q_drv;   // NEW: energy for DRV activation  (J/mol)
  const Real _R_gas_constant; // NEW: gas constant  (such as 8.314 J/mol/K)

  // --- Parameters for Precipitate Strengthening (PDF 3.2) ---
 // g' (L12) parameters
  const Real _C_g_prime;         // NEW: C_gamma' constant in strengthening eqn
  const Real _Gamma_APB_g_prime; // NEW: APB energy for g' (L12)
  const Real _f_vol_g_prime;     // NEW: Volume fraction of g' (L12)

  // g'' (DO22) parameters
  const Real _C_g_pp;            // NEW: C_gamma'' constant in strengthening eqn (pp = prime-prime)
  const Real _Gamma_APB_g_pp;    // NEW: APB energy for g'' (DO22)
  const Real _f_vol_g_pp;        // NEW: Volume fraction of g'' (DO22)

  // --- Initial values for new ISVs ---
  const Real _init_r_eff_g_prime; // NEW: Initial effective radius of g' (L12)
  const Real _init_r_eff_g_pp;    // NEW: Initial effective radius of g'' (DO22)
  const Real _C_shear_g_prime; // Shearing efficiency coefficient for g'
  const Real _C_shear_g_pp;    // Shearing efficiency coefficient for g''

  // Backstress parameters
  const Real _h;
  const Real _h_D;
// --- NEW: Parameters for Intragranular Backstress 
// Agaram, Sukumar, et al. "Dislocation density based crystal plasticity model incorporating the effect of 
// precipitates in IN718 under monotonic and cyclic deformation." 
// International Journal of Plasticity 141 (2021): 102990.
  const Real _k_52; // Controls slope of intra-backstress
  const Real _k_32; // Controls saturation of intra-backstress
  const Real _k_D;  // Controls dislocation generation from precipitates    
  // Initial values of the dislocation density
  const Real _init_rho_ssd;
  const Real _init_rho_gnd_edge;
  const Real _init_rho_gnd_screw;
  
  // Tolerance on dislocation density update
  const Real _rho_tol;
  
  // Dislocation densities
  MaterialProperty<std::vector<Real>> & _rho_ssd;
  const MaterialProperty<std::vector<Real>> & _rho_ssd_old;
  MaterialProperty<std::vector<Real>> & _rho_gnd_edge;
  const MaterialProperty<std::vector<Real>> & _rho_gnd_edge_old;
  MaterialProperty<std::vector<Real>> & _rho_gnd_screw;
  const MaterialProperty<std::vector<Real>> & _rho_gnd_screw_old;

   // --- New ISVs for Precipitate Radius ---
  MaterialProperty<std::vector<Real>> & _r_eff_g_prime; // NEW: Effective radius of g' (ISV)
  const MaterialProperty<std::vector<Real>> & _r_eff_g_prime_old; // NEW
  MaterialProperty<std::vector<Real>> & _r_eff_g_pp;    // NEW: Effective radius of g'' (ISV)
  const MaterialProperty<std::vector<Real>> & _r_eff_g_pp_old;    // NEW

// Backstress variables (SPLIT into Inter- and Intra-)
  MaterialProperty<std::vector<Real>> & _backstress_inter; // RENAMED
  const MaterialProperty<std::vector<Real>> & _backstress_inter_old; // RENAMED
  MaterialProperty<std::vector<Real>> & _backstress_intra; // NEW
  const MaterialProperty<std::vector<Real>> & _backstress_intra_old; // NEW
  
  /// Increment of dislocation densities and backstress
  std::vector<Real> _rho_ssd_increment;
  std::vector<Real> _rho_gnd_edge_increment;
  std::vector<Real> _rho_gnd_screw_increment;
  std::vector<Real> _backstress_inter_increment; // RENAMED
  std::vector<Real> _backstress_intra_increment; // NEW

  std::vector<Real> _r_eff_g_prime_increment; // NEW
  std::vector<Real> _r_eff_g_pp_increment;    // NEW
  /**
   * Stores the values of the dislocation densities and backstress 
   * from the previous substep
   * In classes which use dislocation densities, analogous dislocation density
   * substep vectors will be required.
   */
  std::vector<Real> _previous_substep_rho_ssd;
  std::vector<Real> _previous_substep_rho_gnd_edge;
  std::vector<Real> _previous_substep_rho_gnd_screw;
   std::vector<Real> _previous_substep_backstress_inter; // RENAMED
  std::vector<Real> _previous_substep_backstress_intra; // NEW

  std::vector<Real> _previous_substep_r_eff_g_prime; // NEW
  std::vector<Real> _previous_substep_r_eff_g_pp;    // NEW
  /**
   * Caches the value of the current dislocation densities immediately prior
   * to the update, and they are used to calculate the
   * the dislocation densities for the current substep (or step if
   * only one substep is taken) for the convergence check tolerance comparison.
   */
  std::vector<Real> _rho_ssd_before_update;
  std::vector<Real> _rho_gnd_edge_before_update; 
  std::vector<Real> _rho_gnd_screw_before_update;
  std::vector<Real> _backstress_inter_before_update; // RENAMED
  std::vector<Real> _backstress_intra_before_update; // NEW

  std::vector<Real> _r_eff_g_prime_before_update; // NEW
  std::vector<Real> _r_eff_g_pp_before_update;    // NEW
  /**
   * Flag to include the total twin volume fraction in the plastic velocity
   * gradient calculation, per Kalidindi IJP (2001).
   */
  const bool _include_twinning_in_Lp;

  /**
   * User-defined material property name for the total volume fraction of twins
   * in a twinning propagation constitutive model, when this class is used in
   * conjunction with the twinning propagation model.
   * Note that this value is the OLD material property and thus lags the current
   * value by a single timestep.
   */
  const MaterialProperty<Real> * const _twin_volume_fraction_total;
  
  /**
   * UserObject to read the initial GND density from file
   */
  const ElementPropertyReadFile * const _read_initial_gnd_density;

  /**
   * UserObject to read the grain size for each element from file (New for IN718)
   */
  const ElementPropertyReadFile * const _read_grain_size;  //NEW

  // Directional derivative of the slip rate along the edge dislocation motion direction
  // and along the screw dislocation motion direction
  const bool _include_slip_gradients;
  const ArrayVariableValue & _dslip_increment_dedge;
  const ArrayVariableValue & _dslip_increment_dscrew;
  
  // Temperature dependent properties
  const VariableValue & _temperature;
  const Real _reference_temperature;
  const Real _dCRSS_dT_A;
  const Real _dCRSS_dT_B;
  const Real _dCRSS_dT_C;
  const Real _dCRSS_dT_A_cub; // NEW
  const Real _dCRSS_dT_B_cub; // NEW
  const Real _dCRSS_dT_C_cub;  // NEW
  //Tolerance for gamma_p and gamma_pp radii evolution
  const Real _r_prep_tol;  // NEW 
  // Rotated slip direction to calculate the directional derivative
  // of the plastic strain rate
  // it indicates the edge dislocation velocity direction for all slip systems
  MaterialProperty<std::vector<Real>> & _edge_slip_direction;
  
  // edge dislocation line direction
  // corresponding to direction of motion of screw dislocations
  MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
};
