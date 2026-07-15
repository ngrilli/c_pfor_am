// Nicolò Grilli
// Università di Bristol
// 7 Gennaio 2025

#pragma once

#include "ComputeStressBase.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"

/**
 * A Chaboche model with return mapping
 */
class ChabocheImplicit : public ComputeStressBase
{
public:
  static InputParameters validParams();

  ChabocheImplicit(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpStress();
  
  /// Compute shear and bulk modulus
  virtual void computeElasticConstants();
  
  /// Decompose stress into deviatoric and volumetric
  virtual void decomposeStress(const RankTwoTensor & stress_old);
  
  virtual void decomposeStrainIncrement(const RankTwoTensor & strain_increment);
  
  virtual RankTwoTensor computeTrialStress();
  /**
   * Calculate Mises equivalent stress using the secondInvariant function
   * which is the Mises stress when multiplied by a factor 3
   */
  Real getMisesEquivalent(const RankTwoTensor & stress);
  
  virtual Real yieldFunction(const RankTwoTensor & effective_deviatoric_stress,
                             const Real yield_stress);
                             
  virtual void returnMap(const Real eqvpstrain_old,
                         const RankTwoTensor & plastic_strain_old,
                         Real & eqvpstrain,
                         RankTwoTensor & plastic_strain);
                         
  virtual Real updateIsotropicHardening(const Real eqvpstrain);                 
                         
  /// update backstress after return mapping converged
  virtual RankTwoTensor updateBackstress(const Real delta_gamma,
                                         const RankTwoTensor n,
                                         const unsigned int backstress_index);

  virtual Real numericalJacobian(const Real delta_gamma,
                                 const Real eqvpstrain_old,
                                 const RankTwoTensor n,
                                 const Real residual);
                                
  virtual void elastoPlasticTangentModuli(const Real eqvpstrain,
                                          const RankTwoTensor n);

  // Conversion between symmetric 3x3x3x3 matrix and 6x6 matrix
  virtual std::vector<std::vector<Real>> convertSym3333ToMandel66(const RankFourTensor tensor);
  virtual RankFourTensor convertMandel66ToSym3333(const std::vector<std::vector<Real>> matrix);

  virtual void initMapVoigt();

  virtual void insert66(const std::vector<std::vector<Real>> & A,
                        std::vector<std::vector<Real>> & J,
                        const unsigned int row0,
                        const unsigned int col0);

  // epsilon^p
  MaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
  
  // p
  MaterialProperty<Real> & _eqv_plastic_strain;
  const MaterialProperty<Real> & _eqv_plastic_strain_old;

  // The stress tensor at previous time step
  const MaterialProperty<RankTwoTensor> & _stress_old;

  // delta epsilon_{n+1}
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  
  // Rotation increment in small strain formulation
  const MaterialProperty<RankTwoTensor> & _rotation_increment;
  
  // Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;

  // Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  
  // Backstress tensors
  MaterialProperty<RankTwoTensor> & _backstress1;
  const MaterialProperty<RankTwoTensor> & _backstress1_old;
  MaterialProperty<RankTwoTensor> & _backstress2;
  const MaterialProperty<RankTwoTensor> & _backstress2_old;
  
  // Isotropic hardening R
  MaterialProperty<Real> & _isotropic_hardening;
  const MaterialProperty<Real> & _isotropic_hardening_old;
  
  // Isotropic hardening parameters
  const MaterialProperty<Real> & _sigma_0;
  const MaterialProperty<Real> & _Q;
  const MaterialProperty<Real> & _b;
  
  // Young's modulus and Poisson's ratio as a function of temperature 
  const MaterialProperty<Real> & _E;
  const MaterialProperty<Real> & _nu;
  
  // Backstress parameters
  const MaterialProperty<Real> & _C1;
  const MaterialProperty<Real> & _gamma1;
  const MaterialProperty<Real> & _C2;
  const MaterialProperty<Real> & _gamma2;
  
  // Yield function tolerance
  const Real _tolerance;
  
  // Number of return mapping iterations before unconverged
  const int _max_iterations;

  /// Type of tangent moduli calculation
  const enum class TangentModuliType { ELASTIC, ELASTO_PLASTIC } _tan_mod_type;
  
  // Shear, bulk modulus and first Lamé parameter
  Real _G;
  Real _K;
  Real _lambda;
  
  // Volumetric and deviatoric stress  
  RankTwoTensor _deviatoric_stress_old;
  RankTwoTensor _volumetric_stress_old;
  
  // Volumetric and deviatoric strain increment
  RankTwoTensor _deviatoric_strain_increment;
  RankTwoTensor _volumetric_strain_increment;
  
  // Deviatoric trial stress
  RankTwoTensor _trial_stress;
  
  // Return mapping variables changing at each iteration
  // they represent the quantity at the next time step once the return mapping converged
  RankTwoTensor _effective_deviatoric_stress;
  RankTwoTensor _backstress1_iter;
  RankTwoTensor _backstress2_iter;

  // Mandel-Voigt notation index map and weighting factors
  std::vector<std::vector<unsigned int>> _map_Voigt;
  std::vector<Real> _weight_Mandel;
};
