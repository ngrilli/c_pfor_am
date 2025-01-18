// Nicolò Grilli
// Università di Bristol
// 7 Gennaio 2025

#include "Chaboche.h"

registerMooseObject("c_pfor_amApp", Chaboche);

InputParameters
Chaboche::validParams()
{
  InputParameters params = ComputeStressBase::validParams();

  params.addClassDescription("A Chaboche model with return mapping that follows the implementation in "
                             "0. S. Hopperstad and S. Remseth "
                             "A return mapping algorithm for a class of cyclic plasticity models "
                             "International Journal for Numerical Methods in Engineering, Vol. 38, 549-564 (1995) ");
  params.addRequiredParam<Real>("sigma_0","Constant stress of the isotropic hardening");
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
    _isotropic_hardening_old(getMaterialPropertyOld<Real>(_base_name + "isotropic_hardening"))
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

  // Rotate the stress tensor to the current configuration
  _stress[_qp] = _rotation_increment[_qp] * _stress[_qp] * _rotation_increment[_qp].transpose();

  // Rotate plastic strain tensor to the current configuration
  //_plastic_strain[_qp] =
  //    _rotation_increment[_qp] * _plastic_strain[_qp] * _rotation_increment[_qp].transpose();

  // Calculate the elastic strain_increment
  //_elastic_strain[_qp] = _mechanical_strain[_qp] - _plastic_strain[_qp];

  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}

// Equation (24)
RankTwoTensor
Chaboche::computeTrialStress(const RankTwoTensor & plastic_strain_old,
                             RankTwoTensor & total_strain,
                             const RankFourTensor & E_ijkl)
{
  return E_ijkl * (total_strain - plastic_strain_old);
}
