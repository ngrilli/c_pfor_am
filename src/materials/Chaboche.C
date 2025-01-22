// Nicolò Grilli
// Università di Bristol
// 7 Gennaio 2025

#include "Chaboche.h"
#include "Function.h"

registerMooseObject("c_pfor_amApp", Chaboche);

InputParameters
Chaboche::validParams()
{
  InputParameters params = ComputeStressBase::validParams();

  params.addClassDescription("A Chaboche model with return mapping");
  params.addRequiredParam<Real>("sigma_0","Constant stress of the isotropic hardening");
  params.addRequiredParam<FunctionName>("E","Young's modulus");
  params.addRequiredParam<FunctionName>("nu","Poisson's ratio");
  params.addParam<Real>("tolerance",1e-6,"Yield function tolerance");
  return params;
}

Chaboche::Chaboche(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _plastic_strain(declareProperty<RankTwoTensor>(_base_name + "plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "plastic_strain")),
    _eqv_plastic_strain(declareProperty<Real>(_base_name + "eqv_plastic_strain")),
    _eqv_plastic_strain_old(getMaterialPropertyOld<Real>(_base_name + "eqv_plastic_strain")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _total_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "total_strain")),
    _strain_increment(getMaterialProperty<RankTwoTensor>(_base_name + "strain_increment")),
    _rotation_increment(getMaterialProperty<RankTwoTensor>(_base_name + "rotation_increment")),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _backstress1(declareProperty<RankTwoTensor>(_base_name + "backstress1")),
    _backstress1_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "backstress1")),
    _backstress2(declareProperty<RankTwoTensor>(_base_name + "backstress2")),
    _backstress2_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "backstress2")),
    _isotropic_hardening(declareProperty<Real>(_base_name + "isotropic_hardening")),
    _isotropic_hardening_old(getMaterialPropertyOld<Real>(_base_name + "isotropic_hardening")),
    _E(&this->getFunction("E")),
    _nu(&this->getFunction("nu")),
    _tolerance(getParam<Real>("tolerance"))
{
}

void
Chaboche::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  _plastic_strain[_qp].zero();
  _eqv_plastic_strain[_qp] = 0.0;
}

void
Chaboche::computeQpStress()
{
  // Compute elastic stiffness
  computeElasticConstants();
  
  // Decompose stress 

  // Rotate the stress tensor to the current configuration
  _stress[_qp] = _rotation_increment[_qp] * _stress[_qp] * _rotation_increment[_qp].transpose();

  // Rotate plastic strain tensor to the current configuration
  _plastic_strain[_qp] =
      _rotation_increment[_qp] * _plastic_strain[_qp] * _rotation_increment[_qp].transpose();

  // Calculate the elastic strain_increment
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _plastic_strain[_qp];

  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}

// Compute shear and bulk modulus
void
Chaboche::computeElasticConstants()
{
  _G = _E->value(_t, _q_point[_qp]) / 2.0 / (1.0 + _nu->value(_t, _q_point[_qp]));
  _K = _E->value(_t, _q_point[_qp]) / 3.0 / (1.0 - 2.0 * _nu->value(_t, _q_point[_qp]));
}

// Trial stress assuming all strain increment is elastic
RankTwoTensor
Chaboche::computeTrialStress(const RankTwoTensor & plastic_strain_old,
                             RankTwoTensor & total_strain,
                             const RankFourTensor & E_ijkl)
{
  return E_ijkl * (total_strain - plastic_strain_old);
}

// Calculate Mises equivalent stress
// this calculation intrinsically removes volumetric components
Real
Chaboche::getMisesEquivalent(const RankTwoTensor & stress)
{
  return std::sqrt(3 * stress.secondInvariant());
}

Real
Chaboche::yieldFunction(const RankTwoTensor & stress, 
                        RankTwoTensor & backstress1, 
                        RankTwoTensor & backstress2, 
                        const Real yield_stress)
{
  return getMisesEquivalent(stress - backstress1 - backstress2) - yield_stress;
}
