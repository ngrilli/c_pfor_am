// Nicolò Grilli
// Haziqah Shahari
// University of Bristol
// 25 Marzo 2023

#include "SA508CreepStressUpdate.h"

registerMooseObject("TensorMechanicsApp", SA508CreepStressUpdate);
registerMooseObject("TensorMechanicsApp", ADSA508CreepStressUpdate);

template <bool is_ad>
InputParameters
SA508CreepStressUpdateTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription(
      "This class uses the stress update material in a radial return isotropic power law creep "
      "model. This class can be used in conjunction with other creep and plasticity materials "
      "for more complex simulations. "
      "This model is a modification of the one reported in: "
      "C. Lu et al. "
      "Deformation mechanism-based true-stress creep model for SA508 Gr.3 steel over the temperature range of 450-750C, "
      "Journal of Nuclear Materials 526 (2019) 151776. ");
  params.addCoupledVar("temperature", "Coupled temperature");
  
  // Grain boundary sliding creep rate parameters
  params.addRequiredParam<Real>("A0", "Prefactor of the grain boundary sliding creep rate law");
  params.addRequiredParam<Real>("activation_energy_QA", "Activation energy for grain boundary sliding");
  params.addRequiredParam<Real>("nA_exponent", "Exponent on effective stress in power-law equation for grain boundary sliding");
  params.addParam<Real>("gas_constant", 8.3143, "Universal gas constant");
  
  // Primary creep parameters
  params.addRequiredParam<Real>("beta", "beta material parameter in equations 6a and 6b");
  params.addParam<MaterialPropertyName>("H_name", "H", "H work hardening temperature dependent coefficient");
  
  // Secondary and tertiary creep parameters
  params.addParam<MaterialPropertyName>("k_name", "k", "k temperature dependent coefficient for secondary creep");
  params.addParam<MaterialPropertyName>("Mprime_name", "Mprime", "Mprime temperature dependent tertiary shape parameter");
  
  return params;
}

template <bool is_ad>
SA508CreepStressUpdateTempl<is_ad>::SA508CreepStressUpdateTempl(
    const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _temperature(this->isParamValid("temperature")
                     ? &this->template coupledGenericValue<is_ad>("temperature")
                     : nullptr),               
    _A0(this->template getParam<Real>("A0")),
    _activation_energy_QA(this->template getParam<Real>("activation_energy_QA")),
    _nA_exponent(this->template getParam<Real>("nA_exponent")),
    _gas_constant(this->template getParam<Real>("gas_constant")),
    _beta(this->template getParam<Real>("beta")),
    _H(this->template getMaterialProperty<Real>("H")),
    _k(this->template getMaterialProperty<Real>("k")),
    _Mprime(this->template getMaterialProperty<Real>("Mprime")),
    _exponential(1.0)
{
}

// Assign the value to _exponential for primary creep
template <bool is_ad>
void
SA508CreepStressUpdateTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & /*effective_trial_stress*/,
    const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/)
{
  if (_temperature)
    _exponential = std::exp(-_activation_energy_QA / (_gas_constant * (*_temperature)[_qp]));

}

template <bool is_ad>
template <typename ScalarType>
ScalarType
SA508CreepStressUpdateTempl<is_ad>::computeResidualInternal(
    const GenericReal<is_ad> & effective_trial_stress, const ScalarType & scalar)
{
  // stress_delta is meant to be used as the Mises stress in the creep rate equation 
  const ScalarType stress_delta = effective_trial_stress - _three_shear_modulus * scalar;
  
  // Calculate grain boundary sliding creep rate
  const ScalarType GBS_creep_rate = _A0 * std::pow(stress_delta, _nA_exponent) * _exponential;
  
  // Calculate time scale for primary creep in equation 6c
  ScalarType t_T = 0.0;
  
  if (_beta > 1e-6 && _H[_qp] > 1e-6 && GBS_creep_rate > 1e-12) { // Check that denominator is not too small
	  
    t_T = stress_delta * (_beta - 1.0) / (_beta * _beta * _H[_qp] * GBS_creep_rate);
	  
  } else {
	  
    t_T = 1e12;
	  
  }
  
  // Calculate primary creep rate
  ScalarType primary_creep_rate = 0.0; 
  
  if (_beta > 1.0) {
	
	// Assume that simulation starts at time 0  
    primary_creep_rate = (GBS_creep_rate / (_beta - 1.0)) * std::exp(-_t / t_T);
	  
  }
  
  // Calculate secondary creep rate
  const ScalarType secondary_creep_rate = _k[_qp] * std::exp(_Mprime[_qp] * _k[_qp] * _t);
  
  const ScalarType creep_rate = primary_creep_rate + secondary_creep_rate;

  return creep_rate * _dt - scalar;
}

// Derivative of the creep rate with respect to the stress
// the prefactor _three_shear_modulus and subtraction -1
// is due to the specific formulation in the function
// computeStressDerivative in RadialReturnCreepStressUpdateBase
template <bool is_ad>
GenericReal<is_ad>
SA508CreepStressUpdateTempl<is_ad>::computeDerivative(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - _three_shear_modulus * scalar;
  
  // Calculate grain boundary sliding creep rate
  const GenericReal<is_ad> GBS_creep_rate = _A0 * std::pow(stress_delta, _nA_exponent) * _exponential;
  
  // Derivative of grain boundary sliding creep rate with respect to the stress
  const GenericReal<is_ad> dGBS_creep_rate_dstress = _A0 * _nA_exponent * std::pow(stress_delta, _nA_exponent - 1.0) * _exponential;

  // Calculate the derivative of the time scale for primary creep with respect to the stress
  GenericReal<is_ad> dt_T_dstress = 0.0;

  if (_beta > 1e-6 && _H[_qp] > 1e-6 && GBS_creep_rate > 1e-12) { // Check that denominator is not too small

    dt_T_dstress += (_beta - 1.0) / (_beta * _beta * _H[_qp] * GBS_creep_rate);
  
    dt_T_dstress -= stress_delta * (_beta - 1.0) * dGBS_creep_rate_dstress / (_beta * _beta * _H[_qp] * GBS_creep_rate * GBS_creep_rate);  
  
  }
  
  // Calculate time scale for primary creep in equation 6c
  GenericReal<is_ad> t_T = 0.0;
  
  if (_beta > 1e-6 && _H[_qp] > 1e-6 && GBS_creep_rate > 1e-12) { // Check that denominator is not too small
	  
    t_T = stress_delta * (_beta - 1.0) / (_beta * _beta * _H[_qp] * GBS_creep_rate);
	  
  } else {
	  
    t_T = 1e12;
	  
  }
  
  // Calculate derivative of the creep rate with respect to the stress
  GenericReal<is_ad> creep_rate_derivative = 0.0;
  
  if (_beta > 1.0 && t_T > 1e-6) {
  
    creep_rate_derivative += dGBS_creep_rate_dstress / (_beta - 1.0);
    
    creep_rate_derivative += GBS_creep_rate * _t * dt_T_dstress / ((_beta - 1.0) * t_T * t_T);
    
    creep_rate_derivative *= std::exp(-_t / t_T);
  
  }
  
  // Multiply by -_three_shear_modulus for compatibility with RadialReturnCreepStressUpdateBase
  creep_rate_derivative *= (-1.0 * _three_shear_modulus);
  
  return creep_rate_derivative * _dt - 1.0;
}

template <bool is_ad>
Real
SA508CreepStressUpdateTempl<is_ad>::computeStrainEnergyRateDensity(
    const GenericMaterialProperty<RankTwoTensor, is_ad> & stress,
    const GenericMaterialProperty<RankTwoTensor, is_ad> & strain_rate)
{
  if (_nA_exponent <= 1)
    return 0.0;

  Real creep_factor = _nA_exponent / (_nA_exponent + 1);

  return MetaPhysicL::raw_value(creep_factor * stress[_qp].doubleContraction((strain_rate)[_qp]));
}

template <bool is_ad>
void
SA508CreepStressUpdateTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
SA508CreepStressUpdateTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
SA508CreepStressUpdateTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->_use_substepping != RadialReturnStressUpdateTempl<is_ad>::SubsteppingType::NONE;
}

template class SA508CreepStressUpdateTempl<false>;
template class SA508CreepStressUpdateTempl<true>;
template Real SA508CreepStressUpdateTempl<false>::computeResidualInternal<Real>(const Real &,
                                                                                   const Real &);
template ADReal
SA508CreepStressUpdateTempl<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                      const ADReal &);
template ChainedReal
SA508CreepStressUpdateTempl<false>::computeResidualInternal<ChainedReal>(const Real &,
                                                                            const ChainedReal &);
template ChainedADReal
SA508CreepStressUpdateTempl<true>::computeResidualInternal<ChainedADReal>(const ADReal &,
                                                                             const ChainedADReal &);
