// Zhuohao Song
// Nicolò Grilli
// Università di Bristol
// 26 Maggio 2026

#pragma once

#include "ComputeThermalExpansionEigenstrainBase.h"

/**
 * ComputeInstantaneousMartensiticEigenstrain computes eigenstrain due to martensitic transformation as a function of temperature
 */
template <bool is_ad>
class ComputeInstantaneousMartensiticEigenstrainTempl
  : public ComputeThermalExpansionEigenstrainBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  ComputeInstantaneousMartensiticEigenstrainTempl(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  virtual ValueAndDerivative<is_ad> computeThermalStrain() override;

  const Real _Ms;
  const Real _Mf;
  const Real _martensitic_strain_magnitude;

  /**
   *@{ Stores the martensitic strain as a scalar for use in computing an incremental update.
   */
  GenericMaterialProperty<Real, is_ad> & _martensitic_strain;
  const MaterialProperty<Real> & _martensitic_strain_old;
  //@}

  /// Indicates whether we are on the first step, avoiding false positives when restarting
  bool & _step_one;

  using ComputeThermalExpansionEigenstrainBaseTempl<is_ad>::_qp;
};

typedef ComputeInstantaneousMartensiticEigenstrainTempl<false>
    ComputeInstantaneousMartensiticEigenstrain;
typedef ComputeInstantaneousMartensiticEigenstrainTempl<true>
    ADComputeInstantaneousMartensiticEigenstrain;