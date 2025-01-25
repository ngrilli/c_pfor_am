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
    _total_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "total_strain")), // is this necessary? total or mech?
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
  _backstress1[_qp].zero();
  _backstress2[_qp].zero();
  _isotropic_hardening[_qp] = 0.0; // TO DO: yield surface initialization
}

void
Chaboche::computeQpStress()
{
  // Compute elastic stiffness
  computeElasticConstants();
  
  // Decompose stress
  decomposeStress(_stress_old[_qp]);

  _deviatoric_strain_increment = _strain_increment[_qp].deviatoric();
  
  // Compute trial stress: _trial_stress is deviatoric
  _trial_stress = computeTrialStress(_deviatoric_strain_increment);
  
  _deviatoric_stress.zero();
  _deviatoric_stress = _trial_stress;
  _effective_deviatoric_stress = _deviatoric_stress - _backstress1[_qp] - _backstress2[_qp];
  
  // Return mapping variable that is being solved for
  Real delta_gamma = 0.0;
  
  if (yieldFunction(_effective_deviatoric_stress,_isotropic_hardening[_qp]) > 1.0e-8) {

    // unit vector in stress space that is perpendicular to the yield surface when convergence is reached
    RankTwoTensor n = _effective_deviatoric_stress / getMisesEquivalent(_effective_deviatoric_stress);


    int counter = 0;
    counter++;
  }

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

void
Chaboche::decomposeStress(const RankTwoTensor & stress_old)
{
  _deviatoric_stress_old = stress_old.deviatoric();
  _volumetric_stress_old.zero();
  _volumetric_stress_old.addIa(stress_old.trace() / 3.0);
}

// Trial stress assuming the entire strain increment is elastic
RankTwoTensor
Chaboche::computeTrialStress(const RankTwoTensor & deviatoric_strain_increment)
{
  return _deviatoric_stress_old + 2.0 * _G * deviatoric_strain_increment;
}

// Calculate Mises equivalent stress
// this calculation intrinsically removes volumetric components
Real
Chaboche::getMisesEquivalent(const RankTwoTensor & stress)
{
  return std::sqrt(3 * stress.secondInvariant());
}

Real
Chaboche::yieldFunction(const RankTwoTensor & effective_deviatoric_stress,
                        const Real yield_stress)
{
  return getMisesEquivalent(effective_deviatoric_stress) - yield_stress;
}
