// Nicolo Grilli
// University of Bristol
// Daijun Hu 
// Dai Shi
// National University of Singapore
// 24 Luglio 2021

#pragma once

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticity uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 *
 * Constitutive model from:
 * Nicolo Grilli , Alan C.F. Cocks , Edmund Tarleton
 * Crystal plasticity finite element modelling of coarse-grained alpha-uranium
 * Computational Materials Science 171 (2020) 109276
 *
 */
class FiniteStrainCrystalPlasticityUranium : public FiniteStrainCrystalPlasticity
{
public:
  static InputParameters validParams();

  FiniteStrainCrystalPlasticityUranium(const InputParameters & parameters);

protected:
  /**
   * This function initializes additional parameters.
   * In this case state variables for the dislocation densities
   */
  virtual void initAdditionalProps();
  
  /**
   * This function set variables for internal variable solve.
   * Assign dislocation densities at the previous time step
   * to temporary variables
   */
  virtual void preSolveStatevar();
  
  /**
   * This function update internal variable after solve.
   * Assign updated dislocation density
   * to state variables
   */
  virtual void postSolveStatevar();
  
  /**
   * This function calculates stress residual.
   */
  virtual void calcResidual( RankTwoTensor &resid );
  
  /**
  * This function updates the slip increments.
  * And derivative of slip w.r.t. resolved shear stress.
  */
  virtual void getSlipIncrements();
  
  /**
   * This function updates the slip system resistances
   * based on Taylor hardening model
   */
  virtual void updateGss();
  
  /**
   * This function updates the dislocation densities
   */
  virtual void updateDisloDensity();

  const VariableValue & _temp;
  const Real _thermal_expansion;
  const Real _reference_temperature;

  const Real _dCRSS_dT_A;
  const Real _dCRSS_dT_B;
  const Real _dCRSS_dT_C;
  const Real _dCTE_dT;
  
  const Real _k_bol; // Boltzmann constant, k in equation (8) of the CMS paper
  const Real _da0; // Constant dislocation annihilation length, \hat{d}_\alpha^0 in equation (8)
  
  // Logarithm of the strain rate ratio in equation (8) of the CMS paper
  // \ln \left ( \frac{ \dot{\varepsilon}_0 }{ \dot{\varepsilon} } \right )
  const Real _log_strain_rate_ratio;
  
  const Real _drag_stress; // Drag stress, D_\alpha in equation (8) of the CMS paper  
  
  const Real _ka; // Pre-factor for dislocation multiplication, k_\alpha^1 in equation (6)
  const Real _burgers_vector; // Burgers vector magnitude, b_\alpha in equation (5)
  const Real _projected_mu; // Projected shear modulus on the slip systems, \mu_\alpha in eq. (5)
  const Real _tau0; // Constant friction stress, tau_\alpha^0 in equation (5) of the CMS paper
  
  const Real _init_rho_for; // Initial value of forest dislocation density, same values for all slip systems
  const Real _init_rho_sub; // Initial value of substructure dislocation density
  bool _rho_sub_flag; // Flag to determine whether to include rho_sub in simulations
  
  // critical resolved shear stress
  // exponentially decreased with temperature
  std::vector<Real> _gssT;
  
  // Green-Lagrange strain tensor expressed in the
  // lattice coordinate system
  MaterialProperty<RankTwoTensor> & _lattice_strain;
  
  // Slip increment for output
  MaterialProperty<std::vector<Real>> & _slip_incr_out;
  
  // Forest dislocation density
  MaterialProperty<std::vector<Real>> & _rho_for;
  const MaterialProperty<std::vector<Real>> & _rho_for_old;
  
  // Substructure dislocation density
  MaterialProperty<Real> & _rho_sub;
  const MaterialProperty<Real> & _rho_sub_old;
  
  // Temporary variables for dislocation density
  // dimension of the vectors must be assigned in the constructor
  std::vector<Real> _rho_for_tmp;
  std::vector<Real> _rho_for_tmp_old;
  Real _rho_sub_tmp;
  Real _rho_sub_tmp_old;

};


