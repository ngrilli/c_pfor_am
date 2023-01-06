// Nicolò Grilli
// University of Bristol
// 6 Dicembre 2022

#include "HyperElasticThermalDamage.h"
#include "libmesh/utility.h"

registerMooseObject("TensorMechanicsApp", HyperElasticThermalDamage);

InputParameters
HyperElasticThermalDamage::validParams()
{
  InputParameters params = HyperElasticPhaseFieldIsoDamage::validParams();
  params.addParam<Real>("reference_temperature", 293.0, "Reference temperature at which thermal expansion is zero. ");
  params.addParam<Real>("thermal_expansion_coeff", 0.0, "Thermal expansion coefficient. ");
  params.addRequiredCoupledVar("temperature", "Temperature");
  params.addParam<bool>("suppress_history", false, "Use current elastic energy and do not consider max elastic energy during history. ");
  params.addClassDescription(
      "Computes damaged stress and energy in the intermediate configuration assuming isotropy. "
      "Thermal stress is included. ");

  return params;
}

HyperElasticThermalDamage::HyperElasticThermalDamage(const InputParameters & parameters)
  : HyperElasticPhaseFieldIsoDamage(parameters),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _thermal_expansion_coeff(getParam<Real>("thermal_expansion_coeff")),
    _temperature(coupledValue("temperature")),
    _suppress_history(getParam<bool>("suppress_history"))
{
}

void
HyperElasticThermalDamage::computeDamageStress()
{
  Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);

  Real c = _c[_qp];
  Real xfac = Utility::pow<2>(1.0 - c) + _kdamage;
  
  RankTwoTensor iden(RankTwoTensor::initIdentity);

  std::vector<Real> eigval;
  RankTwoTensor evec;
  _ee.symmetricEigenvaluesEigenvectors(eigval, evec);

  for (const auto i : make_range(Moose::dim))
    _etens[i] = RankTwoTensor::selfOuterProduct(evec.column(i));

  Real etr = 0.0;
  for (const auto i : make_range(Moose::dim))
    etr += eigval[i];

  Real etrpos = (std::abs(etr) + etr) / 2.0;
  Real etrneg = (std::abs(etr) - etr) / 2.0;

  RankTwoTensor pk2pos, pk2neg;

  for (const auto i : make_range(Moose::dim))
  {
    pk2pos += _etens[i] * (lambda * etrpos + 2.0 * mu * (std::abs(eigval[i]) + eigval[i]) / 2.0);
    pk2neg += _etens[i] * (lambda * etrneg + 2.0 * mu * (std::abs(eigval[i]) - eigval[i]) / 2.0);
  }

  _pk2_tmp = pk2pos * xfac - pk2neg;
  
  // Add thermal stress
  _pk2_tmp -= (3*lambda + 2*mu) * _thermal_expansion_coeff * (_temperature[_qp] - _reference_temperature) * iden;

  if (_save_state)
  {
    std::vector<Real> epos(LIBMESH_DIM), eneg(LIBMESH_DIM);
    for (const auto i : make_range(Moose::dim))
    {
      epos[i] = (std::abs(eigval[i]) + eigval[i]) / 2.0;
      eneg[i] = (std::abs(eigval[i]) - eigval[i]) / 2.0;
    }

    // sum squares of epos and eneg
    Real pval(0.0), nval(0.0);
    for (const auto i : make_range(Moose::dim))
    {
      pval += epos[i] * epos[i];
      nval += eneg[i] * eneg[i];
    }

    // Energy with positive principal strains
    const Real G0_pos = lambda * etrpos * etrpos / 2.0 + mu * pval;
    const Real G0_neg = lambda * etrneg * etrneg / 2.0 + mu * nval;

    Real hist_variable = _hist_old[_qp];
    // Assign history variable and derivative
    if (G0_pos > _hist_old[_qp])
      _hist[_qp] = G0_pos;
    else
      _hist[_qp] = _hist_old[_qp];
      
    if (_use_current_hist)
      hist_variable = _hist[_qp];      
      
    // If G0_pos decreases, then decrease the history variable as well
    if (_suppress_history) {
		
      _hist[_qp] = G0_pos;
      hist_variable = G0_pos;
		
	}

    // Elastic free energy density
    _F[_qp] = hist_variable * xfac - G0_neg + _gc[_qp] / (2 * _l[_qp]) * c * c;

    // derivative of elastic free energy density wrt c
    _dFdc[_qp] = -hist_variable * 2.0 * (1.0 - c) * (1 - _kdamage) + _gc[_qp] / _l[_qp] * c;

    // 2nd derivative of elastic free energy density wrt c
    _d2Fdc2[_qp] = hist_variable * 2.0 * (1 - _kdamage) + _gc[_qp] / _l[_qp];

    _dG0_dee = pk2pos;

    _dpk2_dc = -pk2pos * 2.0 * (1.0 - c);
  }
}
