// Nicolo Grilli
// Daijun Hu 
// National University of Singapore
// 16 Novembre 2020

#pragma once

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticityDislo uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Temperature dependence of the CRSS.
 * Dislocation based model.
 * Calculation of the dislocation velocity depending on resolved shear stress.
 * CRSS is calculated using Taylor hardening law and bow-out curvature term for line tension.
 */
class FiniteStrainCrystalPlasticityDislo : public FiniteStrainCrystalPlasticity
{
public:
  static InputParameters validParams();

  FiniteStrainCrystalPlasticityDislo(const InputParameters & parameters);

protected:
  /**
   * This function calculates stress residual.
   */
  virtual void calcResidual( RankTwoTensor &resid );
  
  // Critical resolved shear stress decreases exponentially with temperature
  // A + B exp(- C * (T - 293.0))
  virtual void TempDependCRSS();
  
  /**
  * This function updates the slip increments.
  * And derivative of slip w.r.t. resolved shear stress.
  */
  virtual void getSlipIncrements();
  
  /**
  * This function
  * stores the dislocation velocity value
  * to couple with dislocation transport
  */
  virtual void getDisloVelocity();

  /**
  * This function
  * stores the dislocation velocity direction
  * to couple with dislocation transport
  */
  virtual void OutputSlipDirection();
  
  // Old function: Kept to avoid code break in computeQpStress
  /**
   * This function updates the slip system resistances
   * based on Taylor hardening model
   */
  virtual void updateGss();
  
  /**
   * This function update stress and internal variable after solve.
   */
  virtual void postSolveQp();  
  
  /**
   * This function perform RU decomposition to obtain the rotation tensor.
   */
  RankTwoTensor get_current_rotation(const RankTwoTensor & a);  
  
  /**
   * This function perform RU decomposition to obtain the rotation tensor.
   * Added debug information when RU decomposition fails
   */
  RankTwoTensor getMatRot(const RankTwoTensor & a);

  const VariableValue & _temp;
  
  // coupled curvature (only one slip system, could be extended to 12 by adding more coupled vars)
  const VariableValue & _q_t;
  
  const VariableValue & _rho_edge_pos_1;
  const VariableValue & _rho_edge_neg_1;
  const VariableValue & _rho_edge_pos_2;
  const VariableValue & _rho_edge_neg_2;
  const VariableValue & _rho_edge_pos_3;
  const VariableValue & _rho_edge_neg_3;
  const VariableValue & _rho_edge_pos_4;
  const VariableValue & _rho_edge_neg_4;
  const VariableValue & _rho_edge_pos_5;
  const VariableValue & _rho_edge_neg_5;
  const VariableValue & _rho_edge_pos_6;
  const VariableValue & _rho_edge_neg_6;
  const VariableValue & _rho_edge_pos_7;
  const VariableValue & _rho_edge_neg_7;
  const VariableValue & _rho_edge_pos_8;
  const VariableValue & _rho_edge_neg_8;
  const VariableValue & _rho_edge_pos_9;
  const VariableValue & _rho_edge_neg_9;
  const VariableValue & _rho_edge_pos_10;
  const VariableValue & _rho_edge_neg_10;
  const VariableValue & _rho_edge_pos_11;
  const VariableValue & _rho_edge_neg_11;
  const VariableValue & _rho_edge_pos_12;
  const VariableValue & _rho_edge_neg_12;
  
  const VariableValue & _rho_screw_pos_1;
  const VariableValue & _rho_screw_neg_1;
  const VariableValue & _rho_screw_pos_2;
  const VariableValue & _rho_screw_neg_2;
  const VariableValue & _rho_screw_pos_3;
  const VariableValue & _rho_screw_neg_3;
  const VariableValue & _rho_screw_pos_4;
  const VariableValue & _rho_screw_neg_4;
  const VariableValue & _rho_screw_pos_5;
  const VariableValue & _rho_screw_neg_5;
  const VariableValue & _rho_screw_pos_6;
  const VariableValue & _rho_screw_neg_6;
  const VariableValue & _rho_screw_pos_7;
  const VariableValue & _rho_screw_neg_7;
  const VariableValue & _rho_screw_pos_8;
  const VariableValue & _rho_screw_neg_8;
  const VariableValue & _rho_screw_pos_9;
  const VariableValue & _rho_screw_neg_9;
  const VariableValue & _rho_screw_pos_10;
  const VariableValue & _rho_screw_neg_10;
  const VariableValue & _rho_screw_pos_11;
  const VariableValue & _rho_screw_neg_11;
  const VariableValue & _rho_screw_pos_12;
  const VariableValue & _rho_screw_neg_12;
  
  // Forest dislocation density
  const VariableValue & _rho_forest;
  
  const Real _thermal_expansion;
  const Real _reference_temperature;

  const Real _dCRSS_dT_A;
  const Real _dCRSS_dT_B;
  const Real _dCRSS_dT_C;
  const Real _dislo_mobility;
  const Real _reduced_mobility;
  const Real _burgers_vector_mag;
  const Real _shear_modulus_hardening;
  const Real _dislo_max_velocity;
  
  // Bow-out curvature term for line tension
  // See Hull, Bacon, Dislocations book equation 4.30
  const Real _bowout_coef;
  const Real _bowout_rho_threshold;  
  
  // If the dislocation density goes below the threshold threshold
  // we want the velocity to go to zero because velocity in a region without dislocations
  // is irrelevant for the model and may induce numerical oscillations of the variables
  // This option is activated with the flag _rho_v_thres_flag
  const Real _rho_v_thres;
  bool _rho_v_thres_flag;
  
  // critical resolved shear stress
  // exponentially decreased with temperature
  std::vector<Real> _gssT;
  
  // Rotated slip direction to couple with dislocation transport
  // to indicate dislocation velocity direction for all slip systems
  // edge dislocations
  MaterialProperty<std::vector<Real>> & _edge_slip_direction;
  
  // edge dislocation line direction
  // corresponding to direction of motion of screw dislocations
  MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
  // Slip increment for output
  MaterialProperty<std::vector<Real>> & _slip_incr_out;
  
  // Dislocation velocity
  MaterialProperty<std::vector<Real>> & _dislo_velocity;
  
  // Derivative of the dislocation velocity with respect to the RSS
  // on the same slip system
  MaterialProperty<std::vector<Real>> & _ddislo_velocity_dtau;
  
  // resolved shear stress for output
  // it is used to inform the DoubleCrossSlip object
  MaterialProperty<std::vector<Real>> & _tau_out;

};


