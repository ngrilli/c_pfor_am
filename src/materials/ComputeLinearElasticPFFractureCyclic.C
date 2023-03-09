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
  params.addParam<Real>("tau_cyclic_stress_history", 1.0, "Characteristic time for stress decrease in stress history.");
  params.addParam<MaterialPropertyName>("NS_curve", "NS_curve", "Number of cycles to failure at a given stress level.");
  return params;
}

ComputeLinearElasticPFFractureCyclic::ComputeLinearElasticPFFractureCyclic(
    const InputParameters & parameters)
  : ComputeLinearElasticPFFractureStress(parameters),
  _cycles_per_unit_time(getParam<Real>("cycles_per_unit_time")),
  _tau_cyclic_stress_history(getParam<Real>("tau_cyclic_stress_history")),
  _alpha_cyclic(declareProperty<Real>("alpha_cyclic")),
  _alpha_cyclic_old(getMaterialPropertyOld<Real>("alpha_cyclic")),
  _fatigue_degradation(getMaterialProperty<Real>("fatigue_degradation")),
  _cyclic_stress_history(declareProperty<Real>("cyclic_stress_history")),
  _cyclic_stress_history_old(getMaterialPropertyOld<Real>("cyclic_stress_history")),
  _NS_curve(getMaterialProperty<Real>("NS_curve"))
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
  
  // Eigenvalues of the stress tensor to use in the SN curve
  // and maximum eigenvalue
  RankTwoTensor eigvec_stress;
  std::vector<Real> eigval_stress(LIBMESH_DIM);
  Real max_eigval_stress = 0.0;
  
  // Undamaged stress for fatigue life calculation 
  RankTwoTensor undamaged_stress;

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
  
  // Calculate stress history to determine the point of the SN curve
  // based on the max eigenvalue
  // note that max_eigval_stress cannot become negative
  undamaged_stress = stress0pos - _pressure[_qp] * I2 * _I[_qp] + stress0neg;
  undamaged_stress.positiveProjectionEigenDecomposition(eigval_stress, eigvec_stress);
  
  for (const auto i : make_range(Moose::dim)) {

    max_eigval_stress = std::max(max_eigval_stress,eigval_stress[i]);
	  
  }
  
  // Evolve stress history used for cyclic fatigue calculation
  // _cyclic_stress_history is used in a ParsedMaterial to calculate _NS_curve
  
  if (max_eigval_stress > _cyclic_stress_history_old[_qp]) { // Track maximum stress
	  
    _cyclic_stress_history[_qp] = max_eigval_stress;	  
	  
  } else { // Decrease history variable exponentially in time if max stress is decreasing
	  
    _cyclic_stress_history[_qp] = _cyclic_stress_history_old[_qp] * (1.0 - _dt / _tau_cyclic_stress_history);
	  
  } 

  // Energy with positive principal strains
  F_pos = lambda * etrpos * etrpos / 2.0 + mu * pval;
  F_neg = -lambda * etrneg * etrneg / 2.0 + mu * nval;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];

  _Jacobian_mult[_qp] = (I4sym - (1 - _D[_qp]) * Ppos) * _elasticity_tensor[_qp];

  // update the fatigue effects variable using Miner's rule
  // _NS_curve is used in a ParsedMaterial to calculate _fatigue_degradation
  if (_NS_curve[_qp] == 0)
    mooseError("ComputeLinearElasticPFFractureCyclic: number of cycles to failure must not be zero");
  else
    _alpha_cyclic[_qp] = _alpha_cyclic_old[_qp] + (_cycles_per_unit_time * _dt)/_NS_curve[_qp];

  //fatigue degradation function (1000 its just a value that will be taken from experiment)
  //_fatigue_degradation[_qp] = ((2 * 1000) / (_alpha_cyclic[_qp] + 1000)) * ((2 * 1000) / (_alpha_cyclic[_qp] + 1000));

}

void
ComputeLinearElasticPFFractureCyclic::computeQpStress()
{
  Real F_pos, F_neg;
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  
  // Fpos as modified during cyclic fatigue
  // It is effectively equivalent to a decrease in Gc
  Real cyclic_F_pos;

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
  
  if (_fatigue_degradation[_qp] > 1.0e-6)
    cyclic_F_pos = F_pos / _fatigue_degradation[_qp];
  else 
    cyclic_F_pos = F_pos / 1.0e-6;

  // // Assign history variable
  Real hist_variable = _H_old[_qp];
  if (_use_snes_vi_solver)
  {
    _H[_qp] = cyclic_F_pos;

    if (_use_current_hist)
      hist_variable = _H[_qp];
  }
  else
  {
    if (cyclic_F_pos > _H_old[_qp])
      _H[_qp] = cyclic_F_pos;
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
