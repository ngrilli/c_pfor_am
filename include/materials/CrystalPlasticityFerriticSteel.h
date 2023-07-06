// Michael Salvini
// Nicolò Grilli
// University of Bristol
// 8 Maggio 2023

#pragma once

#include "CrystalPlasticityDislocationUpdateBase.h"
#include "ElementPropertyReadFile.h"

class CrystalPlasticityFerriticSteel;

/**
 * CrystalPlasticityIrradiatedRPVSteel uses the multiplicative decomposition of the
 * deformation gradient and solves the PK2 stress residual equation at the
 * intermediate configuration to evolve the material state. The internal
 * variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Dislocation based model for crystal plasticity with irradiation damage for RPV steel.
 * Additionally GND density due to slip gradients
 * are also included.
 * This model is meant to be used with BCC slip systems only
 */

class CrystalPlasticityFerriticSteel : public CrystalPlasticityDislocationUpdateBase
{
public:
  static InputParameters validParams();

  CrystalPlasticityFerriticSteel(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
   */
  virtual void initQpStatefulProperties() override;

  // Initialize constant reference interaction matrix between slip systems
  virtual void initializeReferenceInteractionMatrix();

  // Initialize the constant slip resistance \tau_0 in equation (3)
  virtual void initializeConstSlipResistance();

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

  virtual bool calculateSlipRate() override;

  virtual void calculateSlipResistance();
  
  virtual void calculateObstaclesDensity();

  virtual void
  calculateEquivalentSlipIncrement(RankTwoTensor & /*equivalent_slip_increment*/) override;



  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;
  
  /*
   * Determines if the state variables, e.g. defect densities, have converged
   * by comparing the change in the values over the iteration period.
   */
  virtual bool areConstitutiveStateVariablesConverged() override;

  /**
   * Stores the current value of the slip system resistance into a separate
   * material property in case substepping is needed.
   */
  virtual void updateSubstepConstitutiveVariableValues() override;

  // Cache the slip system value before the update for the diff in the convergence check
  virtual void cacheStateVariablesBeforeUpdate() override;

  /**
   * Following the Constitutive model for slip system resistance as given in
   * equation (11)
   */
  virtual void calculateStateVariableEvolutionRateComponent() override;

  // Calculate increment of SSD
  virtual void calculateSSDincrement();

  // calculate the irradiation dislocation loops increment based on equation (21)
  virtual void calculateDLincrement();

  // calculate the irradiation solute cluster increment based on equation (23)
  virtual void calculateSCincrement();

  /*
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual bool updateStateVariables() override;

  // model parameters

  // Magnitude of the Burgers vector
  const Real _burgers_vector_mag;

  // Shear modulus in Taylor hardening law G
  // mu in the article
  const Real _shear_modulus;

  // Self and collinear interaction coefficients of the slip systems
  // See figure 1 in
  // Ghiath Monnet, Ludovic Vincent, Lionel Gélébart
  // Multiscale modeling of crystal plasticity in Reactor Pressure Vessel
  // steels: Prediction of irradiation hardening
  // Journal of Nuclear Materials 514 (2019) 128-138
  const Real _a_self;
  const Real _a_col;

  // Constant average diameters of
  // irradiation dislocation loops and
  // irradiation solute clusters
  // and interaction coefficients with the slip systems.
  // According to section 3.2 in
  // Nathan R.Barton, Athanasios Arsenlis, Jaime Marian
  // A polycrystal plasticity model of strain localization in irradiated iron
  // Journal of the Mechanics and Physics of Solids
  // Volume 61, Issue 2, February 2013, Pages 341-351
  // https://www.sciencedirect.com/science/article/pii/S002250961200230X?via%3Dihub
  // the average diameter of the irradiation dislocation loops is 100b
  // we assume the same default value for dislocation loops
  // and solute clusters
  const Real _C_DL_diameter;
  const Real _a_DL;
  const Real _C_SC_diameter;
  const Real _a_SC;

  // prefactor of the irradiation dislocation loops evolution law (adimensional)
  // in equation (21)
  const Real _lambda_DL;

  // prefactor of the irradiation solute cluster evolution law (adimensional)
  // in equation (23)
  const Real _lambda_SC;
  
  // Skip irradiation calculations if false.
  bool _is_irradiated;

  // slip rate coefficient (s^{-1}) for power law slip
  const Real _ao;

  // exponent for slip rate calculated with power law slip
  const Real _xm;
  
  // Activate creep strain rate
  bool _creep_activated;
  
  // Creep rate constants
  const Real _creep_ao;
  const Real _creep_xm;
  const Real _m_exponent;
  const Real _max_stress_ratio;
  const Real _reduced_ao;

  // Constant slip resistances of
  // 110 slip planes
  // 112 slip planes in twinning direction
  // 112 slip planes in anti-twinning direction
  // for the slip systems see Table 1 in:
  // BCC single crystal plasticity modeling and its
  // experimental identification
  // T Yalcinkaya et al 2008 Modelling Simul. Mater. Sci. Eng. 16 085007
  // https://iopscience.iop.org/article/10.1088/0965-0393/16/8/085007/pdf
  //
  // Hall-Petch effect must be included in these constants
  // Carbide-dislocation interaction must be included in these constants
  const Real _const_slip_resistance_110;
  const Real _const_slip_resistance_112_TW;
  const Real _const_slip_resistance_112_AT;
 
  // Dislocation multiplication and annihilation parameters
  const Real _k_0;
  const Real _y_c;

  // Initial values of the dislocation density
  const Real _init_rho_ssd;
  const Real _init_rho_gnd_edge;
  const Real _init_rho_gnd_screw;

  // Initial values of the irradiation defect densities
  const Real _init_C_DL;
  const Real _init_C_SC;

  // Tolerance on dislocation density update
  const Real _rho_tol;

  // State variables
  MaterialProperty<std::vector<Real>> & _rho_ssd;
  const MaterialProperty<std::vector<Real>> & _rho_ssd_old;

  // GND dislocation densities
  MaterialProperty<std::vector<Real>> & _rho_gnd_edge;
  const MaterialProperty<std::vector<Real>> & _rho_gnd_edge_old;
  MaterialProperty<std::vector<Real>> & _rho_gnd_screw;
  const MaterialProperty<std::vector<Real>> & _rho_gnd_screw_old;

  // C_DL: concentration of dislocation loops induced by irradiation
  // on each slip system
  MaterialProperty<std::vector<Real>> & _C_DL;
  const MaterialProperty<std::vector<Real>> & _C_DL_old;

  // C_SC: concentration of solute clusters induced by irradiation
  // on each slip system
  MaterialProperty<std::vector<Real>> & _C_SC;
  const MaterialProperty<std::vector<Real>> & _C_SC_old;

  /// Increment of state variables
  std::vector<Real> _rho_ssd_increment;
  std::vector<Real> _rho_gnd_edge_increment;
  std::vector<Real> _rho_gnd_screw_increment;
  std::vector<Real> _C_DL_increment;
  std::vector<Real> _C_SC_increment;

  /**
   * Stores the values of the state variables from the previous substep
   * In classes which use dislocation densities, analogous dislocation density
   * substep vectors will be required.
   */
  std::vector<Real> _previous_substep_rho_ssd;
  std::vector<Real> _previous_substep_rho_gnd_edge;
  std::vector<Real> _previous_substep_rho_gnd_screw;
  std::vector<Real> _previous_substep_C_DL;
  std::vector<Real> _previous_substep_C_SC;

  /**
   * Caches the value of the current state variables immediately prior
   * to the update, and they are used to calculate the
   * the dislocation densities for the current substep (or step if
   * only one substep is taken) for the convergence check tolerance comparison.
   */
  std::vector<Real> _rho_ssd_before_update;
  std::vector<Real> _rho_gnd_edge_before_update;
  std::vector<Real> _rho_gnd_screw_before_update;
  std::vector<Real> _C_DL_before_update;
  std::vector<Real> _C_SC_before_update;

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

  // Directional derivative of the slip rate along the edge dislocation motion direction
  // and along the screw dislocation motion direction
  const ArrayVariableValue & _dslip_increment_dedge;
  const ArrayVariableValue & _dslip_increment_dscrew;

  // Rotated slip direction to calculate the directional derivative
  // of the plastic strain rate
  // it indicates the edge dislocation velocity direction for all slip systems
  MaterialProperty<std::vector<Real>> & _edge_slip_direction;

  // edge dislocation line direction
  // corresponding to direction of motion of screw dislocations
  MaterialProperty<std::vector<Real>> & _screw_slip_direction;

  // Total density of dislocations including SSD and GND
  // on each slip system
  std::vector<Real> _rho_tot;

  // Total density of local obstacles for each slip system
  std::vector<Real> _rho_obstacles;

  // Constant slip resistance \tau_0
  // for the different slip systems
  // Hall-Petch effect must be included in these constants
  // Carbide-dislocation interaction must be included in these constants
  std::vector<Real> _const_slip_resistance;

  // Reference interaction matrix between slip systems
  // as shown in Figure 1
  DenseMatrix<Real> _a_ref;

};
