// Nicolò Grilli
// University of Bristol
// Daijun Hu 
// Dai Shi
// National University of Singapore
// 10 Luglio 2021

#pragma once

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticity uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 *
 * Constitutive model from:
 * Nicolò Grilli , Alan C.F. Cocks , Edmund Tarleton
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
  * This function
  * stores the dislocation velocity
  * to couple with dislocation transport
  */
  virtual void OutputSlipDirection();
  
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
  
  // critical resolved shear stress
  // exponentially decreased with temperature
  std::vector<Real> _gssT;
  
  // Green-Lagrange strain tensor expressed in the
  // lattice coordinate system
  MaterialProperty<RankTwoTensor> & _lattice_strain; 
  
  // Rotated slip direction to couple with dislocation transport
  // to indicate dislocation velocity direction for all slip systems
  MaterialProperty<std::vector<Real>> & _slip_direction;
  
  // Slip increment for output
  MaterialProperty<std::vector<Real>> & _slip_incr_out;
  
  // Forest dislocation density
  MaterialProperty<std::vector<Real>> & _rho_for;
  const MaterialProperty<std::vector<Real>> & _rho_for_old;
  
  // Substructure dislocation density
  MaterialProperty<Real> & _rho_sub;
  const MaterialProperty<Real> & _rho_sub_old;
  
  // Temporary variables for dislocation density
  std::vector<Real> _rho_for_tmp;
  std::vector<Real> _rho_for_tmp_old;
  Real _rho_sub_tmp;
  Real _rho_sub_tmp_old;

};


