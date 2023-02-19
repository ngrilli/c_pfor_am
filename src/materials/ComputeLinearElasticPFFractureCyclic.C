// Kyprianos Kythreotis
// Nicolò Grilli
// University of Bristol
// 18 Febbraio 2023

#include "ComputeLinearElasticPFFractureCyclic.h"
#include "MathUtils.h"

registerMooseObject("TensorMechanicsApp", ComputeLinearElasticPFFractureCyclic);

InputParameters
ComputeLinearElasticPFFractureCyclic::validParams()
{
  InputParameters params = ComputeLinearElasticPFFractureStress::validParams();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field "
                             "fracture model, with small strain, cyclic model");
  params.addParam<Real>("cycles_per_unit_time", 0, "Number of cycles per unit time.");
  params.addParam<Real>("alpha_cyclic", "Variable accounting for fatigue effects.");
  return params;
}

ComputeLinearElasticPFFractureCyclic::ComputeLinearElasticPFFractureCyclic(
    const InputParameters & parameters)
  : ComputeLinearElasticPFFractureStress(parameters),
  _cycles_per_unit_time(getParam<Real>("cycles_per_unit_time")),
  _alpha_cyclic(declareProperty<Real>("alpha_cyclic")),
  _fatigue_degradation(declareProperty<Real>("fatigue_degradation"))
{
}

void
ComputeLinearElasticPFFractureCyclic::computeStrainSpectral(Real & F_pos, Real & F_neg)
{
  // Isotropic elasticity is assumed and should be enforced
  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);

  RankTwoTensor I2(RankTwoTensor::initIdentity);

  // Compute eigenvectors and eigenvalues of mechanical strain and projection tensor
  RankTwoTensor eigvec;
  std::vector<Real> eigval(LIBMESH_DIM);
  RankFourTensor Ppos =
      _mechanical_strain[_qp].positiveProjectionEigenDecomposition(eigval, eigvec);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  // Calculate tensors of outerproduct of eigen vectors
  std::vector<RankTwoTensor> etens(LIBMESH_DIM);

  for (const auto i : make_range(Moose::dim))
    etens[i] = RankTwoTensor::selfOuterProduct(eigvec.column(i));

  // Separate out positive and negative eigen values
  std::vector<Real> epos(LIBMESH_DIM), eneg(LIBMESH_DIM);
  for (const auto i : make_range(Moose::dim))
  {
    epos[i] = (std::abs(eigval[i]) + eigval[i]) / 2.0;
    eneg[i] = -(std::abs(eigval[i]) - eigval[i]) / 2.0;
  }

  // Seprate positive and negative sums of all eigenvalues
  Real etr = 0.0;
  for (const auto i : make_range(Moose::dim))
    etr += eigval[i];

  const Real etrpos = (std::abs(etr) + etr) / 2.0;
  const Real etrneg = -(std::abs(etr) - etr) / 2.0;

  // Calculate the tensile (postive) and compressive (negative) parts of stress
  RankTwoTensor stress0pos, stress0neg;
  for (const auto i : make_range(Moose::dim))
  {
    stress0pos += etens[i] * (lambda * etrpos + 2.0 * mu * epos[i]);
    stress0neg += etens[i] * (lambda * etrneg + 2.0 * mu * eneg[i]);
  }

  // sum squares of epos and eneg
  Real pval(0.0), nval(0.0);
  for (const auto i : make_range(Moose::dim))
  {
    pval += epos[i] * epos[i];
    nval += eneg[i] * eneg[i];
  }

  _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;

  // Energy with positive principal strains
  F_pos = lambda * etrpos * etrpos / 2.0 + mu * pval;
  F_neg = -lambda * etrneg * etrneg / 2.0 + mu * nval;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];

  _Jacobian_mult[_qp] = (I4sym - (1 - _D[_qp]) * Ppos) * _elasticity_tensor[_qp];

  // update my number of cycles
  _alpha_cyclic[_qp] += _cycles_per_unit_time * _dt;
  //_alpha_cyclic[_qp] = F_pos;

  //fatigue degradation function (1000 its just a value that will be taken from experiment)
  _fatigue_degradation[_qp] = ((2 * 1000) / (_alpha_cyclic[_qp] + 1000)) * ((2 * 1000) / (_alpha_cyclic[_qp] + 1000));

}

void
ComputeLinearElasticPFFractureCyclic::computeQpStress()
{
  Real F_pos, F_neg;
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  switch (_decomposition_type)
  {
    case Decomposition_type::strain_spectral:
      computeStrainSpectral(F_pos, F_neg);
      break;
    case Decomposition_type::strain_vol_dev:
      computeStrainVolDev(F_pos, F_neg);
      break;
    case Decomposition_type::stress_spectral:
      computeStressSpectral(F_pos, F_neg);
      break;
    default:
    {
      RankTwoTensor stress = _elasticity_tensor[_qp] * _mechanical_strain[_qp];
      F_pos = stress.doubleContraction(_mechanical_strain[_qp]) / 2.0;
      F_neg = 0.0;
      if (_use_current_hist)
        _d2Fdcdstrain[_qp] = stress * _dDdc[_qp];

      _stress[_qp] = _D[_qp] * stress - _pressure[_qp] * I2 * _I[_qp];
      _dstress_dc[_qp] = stress * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
      _Jacobian_mult[_qp] = _D[_qp] * _elasticity_tensor[_qp];
    }
  }

  // // Assign history variable
  Real hist_variable = _H_old[_qp];
  if (_use_snes_vi_solver)
  {
    _H[_qp] = F_pos;

    if (_use_current_hist)
      hist_variable = _H[_qp];
  }
  else
  {
    if (F_pos > _H_old[_qp])
      _H[_qp] = F_pos;
    else
      _H[_qp] = _H_old[_qp];

    if (_use_current_hist)
      hist_variable = _H[_qp];

    if (hist_variable < _barrier[_qp])
      hist_variable = _barrier[_qp];
  }

  // Elastic free energy density
  _E[_qp] =
      hist_variable * _D[_qp] + F_neg - _pressure[_qp] * _mechanical_strain[_qp].trace() * _I[_qp];
  _dEdc[_qp] =
      hist_variable * _dDdc[_qp] - _pressure[_qp] * _mechanical_strain[_qp].trace() * _dIdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp] -
                 _pressure[_qp] * _mechanical_strain[_qp].trace() * _d2Id2c[_qp];
}
