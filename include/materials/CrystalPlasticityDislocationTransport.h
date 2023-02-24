// Nicolò Grilli
// University of Bristol
// 25 Gennaio 2023

#pragma once

#include "CrystalPlasticityDislocationUpdateBase.h"
#include "ElementPropertyReadFile.h"

class CrystalPlasticityDislocationTransport;

/**
 * CrystalPlasticityDislocationUpdate uses the multiplicative decomposition of the
 * deformation gradient and solves the PK2 stress residual equation at the
 * intermediate configuration to evolve the material state. The internal
 * variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Dislocation based model for crystal plasticity with Orowan's law.
 * Plastic slip is calculated using the CDD model variables,
 * therefore no update of material properties.
 */

class CrystalPlasticityDislocationTransport : public CrystalPlasticityDislocationUpdateBase
{
public:
  static InputParameters validParams();

  CrystalPlasticityDislocationTransport(const InputParameters & parameters);

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
  
  /**
  * This function calculates the critical resolved shear stress 
  */
  virtual void calculateSlipResistance();
  
  /**
  * This function
  * stores the dislocation velocity value
  * to couple with dislocation transport
  */
  virtual void getDisloVelocity();

  virtual void
  calculateEquivalentSlipIncrement(RankTwoTensor & /*equivalent_slip_increment*/) override;

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  // Cache the slip system value before the update for the diff in the convergence check
  virtual void cacheStateVariablesBeforeUpdate() override;

  virtual void calculateStateVariableEvolutionRateComponent() override;

  /*
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual bool updateStateVariables() override;

  /*
   * Determines if the state variables, e.g. defect densities, have converged
   * by comparing the change in the values over the iteration period.
   * It is a dummy function returning true for this model
   * but necessary because called by ComputeCrystalPlasticityStressDamage
   */
  virtual bool areConstitutiveStateVariablesConverged() override;

  // Temperature dependent properties
  const VariableValue & _temperature;
  
  // Four state variables of the CDD model
  // Each component of the vector variable corresponds
  // to a different slip system
  const VectorVariableValue & _rho_t_vector;
  const VectorVariableValue & _rho_edge_vector;
  const VectorVariableValue & _rho_screw_vector;
  const VectorVariableValue & _q_t_vector;

  // Magnitude of the Burgers vector
  const Real _burgers_vector_mag;
  
  // Shear modulus in Taylor hardening law G
  const Real _shear_modulus;
  
  // Prefactor of Taylor hardening law, alpha
  const Real _alpha_0;
  
  // Latent hardening coefficient
  const Real _r;
  
  // Peierls stress
  const Real _tau_c_0;
  
  // Dislocation mobility and reduced mobility
  // when the velocity goes above the max velocity
  const Real _dislo_mobility;
  const Real _reduced_mobility;
  
  // Max dislocation velocity
  const Real _dislo_max_velocity;
  
  // Bow-out curvature term for line tension
  // See Hull, Bacon, Dislocations book equation 4.30
  const Real _bowout_coef;
  const Real _bowout_rho_threshold;  
  
  // If the dislocation density goes below the threshold
  // we want the velocity to go to zero because velocity in a region without dislocations
  // is irrelevant for the model and may induce numerical oscillations of the variables
  // This option is activated with the flag _rho_v_thres_flag
  const Real _rho_v_thres;
  bool _rho_v_thres_flag;

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
  
  // Temperature dependent critical resolved shear stress
  const Real _reference_temperature;
  const Real _dCRSS_dT_A;
  const Real _dCRSS_dT_B;
  const Real _dCRSS_dT_C;
  
  // Dislocation velocity
  MaterialProperty<std::vector<Real>> & _dislo_velocity;
  
  // Derivative of the dislocation velocity with respect to the RSS
  // on the same slip system
  MaterialProperty<std::vector<Real>> & _ddislo_velocity_dtau;
  
  // Rotated slip direction to calculate the directional derivative
  // of the plastic strain rate
  // it indicates the edge dislocation velocity direction for all slip systems
  MaterialProperty<std::vector<Real>> & _edge_slip_direction;
  
  // edge dislocation line direction
  // corresponding to direction of motion of screw dislocations
  MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
};
