// Nicolò Grilli
// University of Bristol
// 12 Maggio 2022

#pragma once

#include "CrystalPlasticityDislocationUpdateBase.h"
#include "ElementPropertyReadFile.h"

class CrystalPlasticityIrradiatedRPVSteel;

/**
 * CrystalPlasticityIrradiatedRPVSteel uses the multiplicative decomposition of the
 * deformation gradient and solves the PK2 stress residual equation at the
 * intermediate configuration to evolve the material state. The internal
 * variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Dislocation based model for crystal plasticity with irradiation damage for RPV steel.
 * as published in
 * Ghiath Monnet, Ludovic Vincent, Lionel Gelebart
 * Multiscale modeling of crystal plasticity in Reactor Pressure Vessel
 * steels: Prediction of irradiation hardening
 * Journal of Nuclear Materials 514 (2019) 128-138
 * Additionally GND density due to slip gradients
 * are also included
 * This model is meant to be used with BCC slip systems only
 */

class CrystalPlasticityIrradiatedRPVSteel : public CrystalPlasticityDislocationUpdateBase
{
public:
  static InputParameters validParams();

  CrystalPlasticityIrradiatedRPVSteel(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
   */
  virtual void initQpStatefulProperties() override;

  // Initialize constant reference interaction matrix between slip systems
  virtual void initializeReferenceInteractionMatrix();

  // Initialize the constant slip resistance \tau_0 in equation (3)
  virtual void intializeConstSlipResistance();

  // Logarithmic correction to the interaction matrix in equation (7)
  virtual void logarithmicCorrectionInteractionMatrix();

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

  virtual void calculateObstaclesDensity();

  virtual void calculateObstaclesStrength();

  virtual void calculateSelfInteractionSlipResistance();

  virtual void calculateHallPetchSlipResistance();

  virtual void calculateLineTensionSlipResistance();

  virtual bool calculateDragSlipRate();

  virtual bool calculateLatticeFrictionSlipRate();

  virtual void calculateCurvatureDiameter();

  virtual void calculateEffectiveRSS();

  virtual void calculateObstaclesSpacing();

  virtual void calculateAvgLengthScrew();

  virtual void
  calculateEquivalentSlipIncrement(RankTwoTensor & /*equivalent_slip_increment*/) override;

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  virtual void calculateDragSlipRateDerivative();

  virtual void calculateLatticeFrictionSlipRateDerivative();

  // Cache the slip system value before the update for the diff in the convergence check
  virtual void cacheStateVariablesBeforeUpdate() override;

  /**
   * Following the Constitutive model for slip system resistance as given in
   * equation (11)
   */
  virtual void calculateStateVariableEvolutionRateComponent() override;

  // Calculate increment of SSD
  virtual void calculateSSDincrement();

  // calculate dislocation mean free path in equation (19)
  virtual void calculateMeanFreePath();

  // calculate annihilation distance in equation (20)
  virtual void calculateAnnihilationDistance();

  // calculate the irradiation dislocation loops increment based on equation (21)
  virtual void calculateDLincrement();

  // calculate the irradiation solute cluster increment based on equation (23)
  virtual void calculateSCincrement();

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

  // model parameters

  // Magnitude of the Burgers vector
  const Real _burgers_vector_mag;

  // Shear modulus in Taylor hardening law G
  // mu in the article
  const Real _shear_modulus;

  // Shear modulus at room temperature in Taylor hardening law G
  // mu(300 K) in the article
  const Real _RT_shear_modulus;

  // Self and collinear interaction coefficients of the slip systems
  const Real _a_self;
  const Real _a_col;

  // Hall-Petch prefactor
  // N. Tsuchida, H. Masuda, Y. Harada, K. Fukaura, Y. Tomota, K.Nagaid
  // value is reported in reference [22] of the article
  // Effect of ferrite grain size on tensile deformation
  // behavior of a ferrite-cementite low carbon steel
  // Materials Science and Engineering: A
  // Volume 488, Issues 1–2, 15 August 2008, Pages 446-452
  // In figure 7, a value about 480 MPa (micron)^{1/2} is reported
  // at low strain
  const Real _K_Hall_Petch;

  // Average grain size
  const Real _d_grain;

  // Carbide planar density
  // The default value is given by the product C_{carbide} D_{carbide}
  // with values in table 1
  const Real _rho_carbide;

  // Carbide interaction coefficient with the slip systems
  const Real _a_carbide;

  // Constant average diameters of
  // irradiation dislocation loops and
  // irradiation solute clusters
  // and interaction coefficients with the slip systems
  const Real _C_DL_diameter;
  const Real _a_DL;
  const Real _C_SC_diameter;
  const Real _a_SC;

  // According to section 3.2 in
  // Nathan R.Barton, Athanasios Arsenlis, Jaime Marian
  // A polycrystal plasticity model of strain localization in irradiated iron
  // Journal of the Mechanics and Physics of Solids
  // Volume 61, Issue 2, February 2013, Pages 341-351
  // https://www.sciencedirect.com/science/article/pii/S002250961200230X?via%3Dihub
  // the average diameter of the irradiation dislocation loops is 100b
  // we assume the same default value for dislocation loops
  // and solute clusters

  // Reference dislocation density at which the interaction
  // matrix between slip system is the reference matrix
  // (1 / micron^2)
  const Real _rho_ref;

  // slip rate coefficient (s^{-1}) in equation (2)
  const Real _ao;

  // exponent for slip rate in equation (2)
  const Real _xm;

  // attack frequency for the lattice friction slip rate (1/s)
  const Real _attack_frequency;

  // minimum length of screw dislocation segments l_c in equation (10)
  const Real _minimum_screw_length;

  // Gibbs free energy jump for thermal slip in equation (3)
  const Real _Gibbs_free_energy_slip;

  // Boltzmann constant
  const Real _k;

  // Constant slip resistances of
  // 110 slip planes
  // 112 slip planes in twinning direction
  // 112 slip planes in anti-twinning direction
  // for the slip systems see Table 1 in:
  // BCC single crystal plasticity modeling and its
  // experimental identification
  // T Yalcinkaya et al 2008 Modelling Simul. Mater. Sci. Eng. 16 085007
  // https://iopscience.iop.org/article/10.1088/0965-0393/16/8/085007/pdf
  const Real _const_slip_resistance_110;
  const Real _const_slip_resistance_112_TW;
  const Real _const_slip_resistance_112_AT;

  // number of intersections with primary and forest dislocations before immobilization
  // in equation (19)
  const Real _K_self;
  const Real _K_forest;

  // annihilation distance that prevails at high temperature in the drag regime
  // in equation (20)
  const Real _y_drag;

  // prefactor of the irradiation dislocation loops evolution law (adimensional)
  // in equation (21)
  const Real _lambda_DL;

  // prefactor of the irradiation solute cluster evolution law (adimensional)
  // in equation (23)
  const Real _lambda_SC;

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

  // _rho_ssd corresponds to rho_s in equation (18)
  MaterialProperty<std::vector<Real>> & _rho_ssd;
  const MaterialProperty<std::vector<Real>> & _rho_ssd_old;

  // GND dislocation densities: not in the original model
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

  // Temperature in K as coupled variables
  // so temperature evolution can be included in the model
  const VariableValue & _temperature;
  
  // Use the Gibbs energy based term for slip
  // If false, a simple power law slip equation is used
  bool _use_lattice_friction_slip;

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

  // Average obstacles strength for each slip system
  // according to equation (6)
  std::vector<Real> _obstacles_strength;

  // Self interaction stress tau_self for each slip system
  std::vector<Real> _tau_self;

  // Hall-Petch slip resistance
  // it is the same for each slip system
  // note that if GND are activated, part of the Hall-Petch effect
  // will be provided by the GNDs
  Real _tau_Hall_Petch;

  // Line tension slip resistance
  std::vector<Real> _tau_line_tension;

  // Drag contribution to the slip increment
  // according to equation (2)
  std::vector<Real> _drag_slip_increment;

  // and its derivative with respect to the RSS
  std::vector<Real> _ddrag_slip_increment_dtau;

  // Lattice friction contribution to the slip increment
  // according to equation (3)
  std::vector<Real> _lattice_friction_slip_increment;

  // and its derivative with respect to the RSS
  std::vector<Real> _dlattice_friction_slip_increment_dtau;

  // Lambda^s in equation (18)
  std::vector<Real> _dislocation_mean_free_path;

  // Slip system depedent annihilation distance
  // in equation (20)
  std::vector<Real> _annihilation_distance;

  // Low mobility of screw dislocations induces a curvature of non-screw
  // dislocations given by the diameter in equation (8)
  std::vector<Real> _curvature_diameter;

  // effective resolved shear stress \tau_{eff}^s
  // for each slip system in equation (4)
  std::vector<Real> _effective_RSS;

  // obstacles spacing in equation (9)
  // \lambda^s
  std::vector<Real> _obstacles_spacing;

  // average length of screw dislocations l_{sc}^s
  // in equation (10)
  std::vector<Real> _avg_length_screw;

  // Constant slip resistance \tau_0
  // for the different slip systems in equation (3)
  std::vector<Real> _const_slip_resistance;

  // Reference interaction matrix between slip systems
  // as shown in Figure 1
  DenseMatrix<Real> _a_ref;

  // Corrected interaction matrix between slip systems
  // that accounts for the logarithmic correction in equation (7)
  DenseMatrix<Real> _a_slip_slip_interaction;

};
