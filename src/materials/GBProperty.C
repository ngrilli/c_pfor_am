// Nicolò Grilli
// Università di Bristol
// 21 Agosto 2024

#include "GBProperty.h"

template <bool is_ad>
InputParameters
GBPropertyTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the isotropic grain growth model "
      "using the same convention as "
      "Min Yang, Lu Wang, Wentao Yan "
      "Phase-ﬁeld modeling of grain evolutions in additive manufacturing from nucleation, growth, to coarsening "
      "npj Computational Materials (2021) 7:56. ");
  params.addRequiredCoupledVar("T", "Temperature in Kelvin");
  params.addParam<Real>("Delta_f_g", 0.125, "Maximum height of the barrier in the free energy density for grains. ");
  params.addParam<Real>("sigma_g_0", 0.385, "Isotropic grain boundary energy in J/m^2. ");
  params.addParam<Real>("Q_g", 140.0, "Grain boundary migration activation energy in kJ/mol. ");
  params.addParam<Real>("D0", 0.034, "Pre-exponential coefﬁcient for the grain boundary mobility prefactor in m^4/(J*s). ");
  params.addParam<Real>("l_g", 9.6, "Grain boundary width in micron. ");
  return params;
}

template <bool is_ad>
GBPropertyTempl<is_ad>::GBPropertyTempl(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _T(coupledValue("T")),
    _Delta_f_g(getParam<Real>("Delta_f_g")),
    _sigma_g_0(getParam<Real>("sigma_g_0")),
    _Q_g(getParam<Real>("Q_g")),
    _D0(getParam<Real>("D0")),
    _l_g(getParam<Real>("l_g")),

    // sigma_g in equation (12), isotropic part
    _sigma(declareGenericProperty<Real, is_ad>("sigma")),
    
    // D_g in equation (18)
    _M_GB(declareGenericProperty<Real, is_ad>("M_GB")),
    
    // k_g in equation (16)
    _kappa(declareGenericProperty<Real, is_ad>("kappa_op")),
    
    // gamma parameter in equation (10)
    _gamma(declareGenericProperty<Real, is_ad>("gamma_asymm")),
    
    // L_g in equation (17) and its derivative
    // with respect to temperature
    _L(declareGenericProperty<Real, is_ad>("L")),
    _dLdT(parameters.hasDefaultCoupledValue("T")
              ? nullptr
              : &declarePropertyDerivative<Real>("L", coupledName("T", 0))),
              
    // l_g in equations (14) and (16)
    _l_GB(declareGenericProperty<Real, is_ad>("l_GB")),
    
    // m_g in equation (14)
    _mu(declareGenericProperty<Real, is_ad>("mu")),
    
    // unused properties
    _entropy_diff(declareGenericProperty<Real, is_ad>("entropy_diff")),
    _molar_volume(declareGenericProperty<Real, is_ad>("molar_volume")),
    _act_wGB(declareGenericProperty<Real, is_ad>("act_wGB")),
    
    // Boltzmann constant in kJ/mol/K
    _kb(8.314462618e-3)
{
}

template <bool is_ad>
void
GBPropertyTempl<is_ad>::computeQpProperties()
{
  // GB mobility derivative
  Real dM_GBdT;
  
  // D_g in equation (18)
  _M_GB[_qp] = _D0 * std::exp(-_Q_g / (_kb * _T[_qp]));
  dM_GBdT = MetaPhysicL::raw_value(_M_GB[_qp] * _Q_g / (_kb * _T[_qp] * _T[_qp]));

  // l_g in equations (14) and (16)
  _l_GB[_qp] = _l_g;

  // L_g in equation (17)
  _L[_qp] = 4.0 / 3.0 * _M_GB[_qp] / _l_g;
  if (_dLdT)
    (*_dLdT)[_qp] = MetaPhysicL::raw_value(4.0 / 3.0 * dM_GBdT / _l_g);
    
  // sigma_g in equation (12), isotropic part
  _sigma[_qp] = _sigma_g_0;

  // equation (16) for k_g
  _kappa[_qp] = 3.0 / 4.0 * _sigma[_qp] * _l_g;
  
  // gamma parameter in equation (10)
  _gamma[_qp] = 1.5;
  
  // equation (14) for m_g
  _mu[_qp] = 3.0 / 4.0 * 1.0 / _Delta_f_g * _sigma_g_0 / _l_g;
  
  // unused properties
  _entropy_diff[_qp] = 0.0;
  _molar_volume[_qp] = 0.0;
  _act_wGB[_qp] = 0.0;
}

template class GBPropertyTempl<false>;
template class GBPropertyTempl<true>;
