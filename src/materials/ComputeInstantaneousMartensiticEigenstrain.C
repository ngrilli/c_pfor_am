// Zhuohao Song
// Nicolò Grilli
// Università di Bristol
// 26 Maggio 2026

#include "ComputeInstantaneousMartensiticEigenstrain.h"
#include "CastDualNumber.h"

registerMooseObject("c_pfor_amApp", ComputeInstantaneousMartensiticEigenstrain);
registerMooseObject("c_pfor_amApp", ADComputeInstantaneousMartensiticEigenstrain);

template <bool is_ad>
InputParameters
ComputeInstantaneousMartensiticEigenstrainTempl<is_ad>::validParams()
{
  InputParameters params = ComputeThermalExpansionEigenstrainBaseTempl<is_ad>::validParams();
  params.addClassDescription("Computes eigenstrain due to martensitic transformation as a function of temperature");
  params.addParam<Real>("Ms", 10000.0, "Martensitic start temperature");
  params.addParam<Real>("Mf", 0.0, "Martensitic finish temperature");
  params.addParam<Real>("martensitic_strain_magnitude", 0.0, "Magnitude of the linear martensitic transformation strain");
  params.suppressParameter<bool>("use_old_temperature");
  return params;
}

template <bool is_ad>
ComputeInstantaneousMartensiticEigenstrainTempl<is_ad>::ComputeInstantaneousMartensiticEigenstrainTempl(const InputParameters & parameters)
  : ComputeThermalExpansionEigenstrainBaseTempl<is_ad>(parameters),
    _Ms(this->template getParam<Real>("Ms")),
    _Mf(this->template getParam<Real>("Mf")),
    _martensitic_strain_magnitude(this->template getParam<Real>("martensitic_strain_magnitude")),
    _martensitic_strain(this->template declareGenericProperty<Real, is_ad>(
        this->_base_name + "InstantaneousMartensitic_strain")),
    _martensitic_strain_old(this->template getMaterialPropertyOld<Real>(
        this->_base_name + "InstantaneousMartensitic_strain")),
    _step_one(this->template declareRestartableData<bool>("step_one", true))
{
  if (this->_use_old_temperature)
    this->paramError("use_old_temperature",
                     "The old temperature value cannot be used in this incremental update model.");
}

template <bool is_ad>
void
ComputeInstantaneousMartensiticEigenstrainTempl<is_ad>::initQpStatefulProperties()
{
  _martensitic_strain[_qp] = 0;
}

template <bool is_ad>
ValueAndDerivative<is_ad>
ComputeInstantaneousMartensiticEigenstrainTempl<is_ad>::computeThermalStrain()
{
  if (this->_t_step > 1)
    _step_one = false;

  const auto & old_temp =
      (_step_one ? this->_stress_free_temperature[_qp] : this->_temperature_old[_qp]);
  const auto delta_T = this->_temperature[_qp] - old_temp;
  
  auto martensitic_strain = _martensitic_strain_old[_qp] + 0.0 * delta_T; // Initialize as ADReal type

  if (delta_T < 0.0) { // cooling stage: martensitic strain is expansive
    martensitic_strain = _martensitic_strain_old[_qp] - _martensitic_strain_magnitude * delta_T / (_Ms - _Mf);
  }

  // Store martensitic strain for use in the next timestep (casts ValueAndDerivative<is_ad>
  // to GenericReal<is_ad>).
  _martensitic_strain[_qp] = dual_number_cast<GenericReal<is_ad>>(martensitic_strain);

  return martensitic_strain;
}

template class ComputeInstantaneousMartensiticEigenstrainTempl<false>;
template class ComputeInstantaneousMartensiticEigenstrainTempl<true>;
