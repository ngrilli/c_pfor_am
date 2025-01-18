// Nicolò Grilli
// Università di Bristol
// 7 Gennaio 2025

#pragma once

#include "ComputeStressBase.h"

/**
 * A Chaboche model with return mapping that follows the implementation in
 * 0. S. Hopperstad and S. Remseth
 * A return mapping algorithm for a class of cyclic plasticity models
 * International Journal for Numerical Methods in Engineering, Vol. 38, 549-564 (1995)
 */
class Chaboche : public ComputeStressBase
{
public:
  static InputParameters validParams();

  Chaboche(const InputParameters & parameters);

protected:
  virtual void computeQpStress();
  virtual void initQpStatefulProperties();
  
  virtual RankTwoTensor computeTrialStress(const RankTwoTensor & plastic_strain_old,
                                           RankTwoTensor & total_strain,
                                           const RankFourTensor & E_ijkl);

  // epsilon^p
  MaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
  
  // p
  MaterialProperty<Real> & _eqv_plastic_strain;
  const MaterialProperty<Real> & _eqv_plastic_strain_old;

  // The stress tensor at previous time step
  const MaterialProperty<RankTwoTensor> & _stress_old;
  
  // epsilon_n
  const MaterialProperty<RankTwoTensor> & _total_strain_old;

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
};
