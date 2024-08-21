// Nicolò Grilli
// Università di Bristol
// 21 Agosto 2024

#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

/// Computes necessary material properties for the isotropic grain growth model "
/// using the same convention as "
/// Min Yang, Lu Wang, Wentao Yan "
/// Phase-ﬁeld modeling of grain evolutions in additive manufacturing from nucleation, growth, to coarsening "
/// npj Computational Materials (2021) 7:56. "

template <bool is_ad>
class GBPropertyTempl : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  GBPropertyTempl(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const VariableValue & _T;

  const Real _Delta_f_g;
  const Real _sigma_g_0;
  const Real _Q_g;
  const Real _D0;
  const Real _l_g;

  // sigma_g in equation (12), isotropic part
  GenericMaterialProperty<Real, is_ad> & _sigma;
  
  // D_g in equation (18)
  GenericMaterialProperty<Real, is_ad> & _M_GB;
  
  // equation (16) for k_g
  GenericMaterialProperty<Real, is_ad> & _kappa;
  
  // gamma parameter in equation (10)
  GenericMaterialProperty<Real, is_ad> & _gamma;
  
  // L_g in equation (17) and its derivative
  // with respect to temperature
  GenericMaterialProperty<Real, is_ad> & _L;
  MaterialProperty<Real> * _dLdT;
  
  // l_g in equations (14) and (16)
  GenericMaterialProperty<Real, is_ad> & _l_GB;
  
  // m_g in equation (14)
  GenericMaterialProperty<Real, is_ad> & _mu;
  
  // unused properties
  GenericMaterialProperty<Real, is_ad> & _entropy_diff;
  GenericMaterialProperty<Real, is_ad> & _molar_volume;
  GenericMaterialProperty<Real, is_ad> & _act_wGB;

  // Boltzmann constant in kJ/mol/K
  const Real _kb;
};

typedef GBPropertyTempl<false> GBProperty;
typedef GBPropertyTempl<true> ADGBProperty;
