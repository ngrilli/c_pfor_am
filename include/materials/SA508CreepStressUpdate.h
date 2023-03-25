// Nicolò Grilli
// Haziqah Shahari
// University of Bristol
// 25 Marzo 2023

#pragma once

#include "RadialReturnCreepStressUpdateBase.h"

/**
 * This class uses the stress update material in a radial return isotropic creep
 * model.  This class is one of the basic radial return constitutive models; more complex
 * constitutive models combine creep and plasticity.
 *
 * This class inherits from RadialReturnCreepStressUpdateBase and must be used
 * in conjunction with ComputeMultipleInelasticStress.  This class calculates
 * creep based on stress, temperature, and time effects.  This class also
 * computes the creep strain as a stateful material property.
 * 
 * This model is a modification of the one reported in:
 * C. Lu et al.
 * Deformation mechanism-based true-stress creep model for SA508 Gr.3 steel over the temperature range of 450-750C,
 * Journal of Nuclear Materials 526 (2019) 151776.
 */
template <bool is_ad>
class SA508CreepStressUpdateTempl : public RadialReturnCreepStressUpdateBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  SA508CreepStressUpdateTempl(const InputParameters & parameters);

  virtual Real computeStrainEnergyRateDensity(
      const GenericMaterialProperty<RankTwoTensor, is_ad> & stress,
      const GenericMaterialProperty<RankTwoTensor, is_ad> & strain_rate) override;

  virtual bool substeppingCapabilityEnabled() override;

  virtual void resetIncrementalMaterialProperties() override;

protected:
  virtual void
  computeStressInitialize(const GenericReal<is_ad> & effective_trial_stress,
                          const GenericRankFourTensor<is_ad> & elasticity_tensor) override;

  virtual void
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & plastic_strain_increment) override;

  virtual GenericReal<is_ad> computeResidual(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar) override
  {
    return computeResidualInternal<GenericReal<is_ad>>(effective_trial_stress, scalar);
  }
  virtual GenericReal<is_ad> computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                               const GenericReal<is_ad> & scalar) override;
  virtual GenericChainedReal<is_ad>
  computeResidualAndDerivative(const GenericReal<is_ad> & effective_trial_stress,
                               const GenericChainedReal<is_ad> & scalar) override
  {
    return computeResidualInternal<GenericChainedReal<is_ad>>(effective_trial_stress, scalar);
  }

  /// Temperature variable value
  const GenericVariableValue<is_ad> * const _temperature;

  /// Prefactor of the grain boundary sliding creep rate law
  const Real _A0;
  
  /// Activation energy for grain boundary sliding
  const Real _activation_energy_QA;
  
  /// Exponent on effective stress in power-law equation for grain boundary sliding
  const Real _nA_exponent;
  
  /// Gas constant for exp term
  const Real _gas_constant;

  /// beta material parameter in equations 6a and 6b
  const Real _beta;
  
  /// H work hardening temperature dependent coefficient
  const MaterialProperty<Real> & _H;
  
  /// k temperature dependent coefficient for secondary creep
  /// see equation 6a
  const MaterialProperty<Real> & _k;
  
  /// Mprime tertiary shape parameter
  /// see equation 6a
  const MaterialProperty<Real> & _Mprime;

  /// Exponential calculated from activation, gas constant, and temperature
  GenericReal<is_ad> _exponential;

  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_qp;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_dt;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_t;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_three_shear_modulus;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_creep_strain;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_creep_strain_old;

private:
  template <typename ScalarType>
  ScalarType computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                     const ScalarType & scalar);
};

typedef SA508CreepStressUpdateTempl<false> SA508CreepStressUpdate;
typedef SA508CreepStressUpdateTempl<true> ADSA508CreepStressUpdate;
