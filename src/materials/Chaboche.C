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
  params.addRequiredParam<FunctionName>("sigma_0","Constant stress of the isotropic hardening");
  params.addRequiredParam<FunctionName>("Q","Maximum isotropic hardening");
  params.addRequiredParam<FunctionName>("b","Hardening rate");
  params.addRequiredParam<FunctionName>("E","Young's modulus");
  params.addRequiredParam<FunctionName>("nu","Poisson's ratio");
  params.addRequiredParam<FunctionName>("C1","C constant for the first backstress");
  params.addRequiredParam<FunctionName>("gamma1","gamma constant for the first backstress");
  params.addRequiredParam<FunctionName>("C2","C constant for the second backstress");
  params.addRequiredParam<FunctionName>("gamma2","gamma constant for the second backstress");
  params.addParam<Real>("tolerance",1e-6,"Yield function tolerance");
  params.addParam<int>("max_iterations",1000,"Number of return mapping iterations before unconverged");
  return params;
}

Chaboche::Chaboche(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _plastic_strain(declareProperty<RankTwoTensor>(_base_name + "plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "plastic_strain")),
    _eqv_plastic_strain(declareProperty<Real>(_base_name + "eqv_plastic_strain")),
    _eqv_plastic_strain_old(getMaterialPropertyOld<Real>(_base_name + "eqv_plastic_strain")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
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
    _sigma_0(&this->getFunction("sigma_0")),
    _Q(&this->getFunction("Q")),
    _b(&this->getFunction("b")),
    _E(&this->getFunction("E")),
    _nu(&this->getFunction("nu")),
    _C1(&this->getFunction("C1")),
    _gamma1(&this->getFunction("gamma1")),
    _C2(&this->getFunction("C2")),
    _gamma2(&this->getFunction("gamma2")),
    _tolerance(getParam<Real>("tolerance")),
    _max_iterations(getParam<int>("max_iterations"))
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
  _isotropic_hardening[_qp] = _sigma_0->value(_t, _q_point[_qp]);
}

void
Chaboche::computeQpStress()
{
  // Compute elastic stiffness
  computeElasticConstants();
  
  // Decompose stress
  decomposeStress(_stress_old[_qp]);
  
  // Decompose strain increment
  decomposeStrainIncrement(_strain_increment[_qp]);

  // Compute trial stress: _trial_stress is deviatoric
  _trial_stress = computeTrialStress();
  
  // Assume first this strain increment does not induce any plasticity
  _eqv_plastic_strain[_qp] = _eqv_plastic_strain_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
  _backstress1[_qp] = _backstress1_old[_qp];
  _backstress2[_qp] = _backstress2_old[_qp];
  _isotropic_hardening[_qp] = _isotropic_hardening_old[_qp];  
  
  // Compute effective stress
  _effective_deviatoric_stress = _trial_stress - _backstress1[_qp] - _backstress2[_qp]; 
  
  if (yieldFunction(_effective_deviatoric_stress,_isotropic_hardening[_qp]) >= 0.0) {
    
    returnMap(_eqv_plastic_strain_old[_qp],
              _plastic_strain_old[_qp],
              _eqv_plastic_strain[_qp],
              _plastic_strain[_qp]);
    
  } else { // purely elastic stress update

    _stress[_qp] = _stress_old[_qp] + 2.0 * _G * _volumetric_strain_increment;
    _stress[_qp].addIa(_lambda * _volumetric_strain_increment.trace());
    _stress[_qp] += 2.0 * _G * _deviatoric_strain_increment;
  }

  // Rotate the stress tensor to the current configuration
  _stress[_qp] = _rotation_increment[_qp] * _stress[_qp] * _rotation_increment[_qp].transpose();

  // Rotate plastic strain tensor to the current configuration
  _plastic_strain[_qp] =
      _rotation_increment[_qp] * _plastic_strain[_qp] * _rotation_increment[_qp].transpose();

  // Calculate the elastic strain
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _plastic_strain[_qp];

  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}

// Compute shear and bulk modulus
void
Chaboche::computeElasticConstants()
{
  _G = _E->value(_t, _q_point[_qp]) / 2.0 / (1.0 + _nu->value(_t, _q_point[_qp]));
  _K = _E->value(_t, _q_point[_qp]) / 3.0 / (1.0 - 2.0 * _nu->value(_t, _q_point[_qp]));
  _lambda = _K - (2.0 / 3.0) * _G;
}

void
Chaboche::decomposeStress(const RankTwoTensor & stress_old)
{
  _deviatoric_stress_old = stress_old.deviatoric();
  _volumetric_stress_old.zero();
  _volumetric_stress_old.addIa(stress_old.trace() / 3.0);
}

void
Chaboche::decomposeStrainIncrement(const RankTwoTensor & strain_increment)
{
  _deviatoric_strain_increment = strain_increment.deviatoric();
  _volumetric_strain_increment.zero();
  _volumetric_strain_increment.addIa(strain_increment.trace() / 3.0);
}

// Trial stress assuming the entire strain increment is elastic
// Trial stress includes only the deviatoric component
RankTwoTensor
Chaboche::computeTrialStress()
{
  return _deviatoric_stress_old + 2.0 * _G * _deviatoric_strain_increment;
}

// Calculate Mises equivalent stress
// this calculation intrinsically removes volumetric components
Real
Chaboche::getMisesEquivalent(const RankTwoTensor & stress)
{
  return std::sqrt(3 * stress.secondInvariant());
}

// Yield function, plasticity if positive
Real
Chaboche::yieldFunction(const RankTwoTensor & effective_deviatoric_stress,
                        const Real yield_stress)
{
  return getMisesEquivalent(effective_deviatoric_stress) - yield_stress;
}

void
Chaboche::returnMap(const Real eqvpstrain_old,
                    const RankTwoTensor & plastic_strain_old,
                    Real & eqvpstrain,
                    RankTwoTensor & plastic_strain)
{
  // scalar plastic multiplier that is solved for
  Real delta_gamma = 0.0;
  
  // Residual and Jacobian as a function of delta_gamma
  Real residual;
  Real jacobian;
  
  // Direction of plastic flow, in tensor form
  RankTwoTensor n;
  
  // Scalar plastic strain increment
  Real delta_eps_p;
  
  // Mises equivalent of the effective deviatoric stress
  Real eqv_effective_deviatoric_stress;
  
  // Isotropic hardening parameters
  Real Q = _Q->value(_t, _q_point[_qp]);
  Real b = _b->value(_t, _q_point[_qp]);
  
  // Iteration counter
  unsigned int i = 0;
  
  do {
    // update direction of plastic flow 
    eqv_effective_deviatoric_stress = getMisesEquivalent(_effective_deviatoric_stress);
    n = _effective_deviatoric_stress / eqv_effective_deviatoric_stress;
    
    // Update plastic strain increment
    delta_eps_p = std::sqrt(2.0 / 3.0) * std::abs(delta_gamma);    
    
    // Update isotropic hardening
    eqvpstrain = eqvpstrain_old + delta_eps_p;
    updateIsotropicHardening(eqvpstrain);
    
    // Calculate residual and Jacobian
    residual = eqv_effective_deviatoric_stress - 3 * _G * delta_gamma - _isotropic_hardening[_qp];
    jacobian = -3 * _G - Q * b * std::exp(-b * eqvpstrain) * std::sqrt(2.0 / 3.0);
    
    // Update plastic multiplier
    delta_gamma -= residual / jacobian;
    
    // Update stress variables
    _effective_deviatoric_stress = _trial_stress - _backstress1[_qp] - _backstress2[_qp] - 2.0 * _G * delta_gamma * n;

    if (i > _max_iterations) // unconverged
      mooseError("Constitutive failure");
    
    // update iteration counter
    i++;
  }
  while (residual > _tolerance);
  
  // Update plastic strain, isotropic hardening and back stress
  delta_eps_p = std::sqrt(2.0 / 3.0) * std::abs(delta_gamma);
  eqvpstrain = eqvpstrain_old + delta_eps_p;
  updateIsotropicHardening(eqvpstrain);
  plastic_strain = plastic_strain_old + delta_gamma * n;
  updateBackstress(delta_gamma,n);
  
  // Update total stress: recombine deviatoric and hydrostatic parts
  _stress[_qp] = _stress_old[_qp] + 2.0 * _G * _volumetric_strain_increment;
  _stress[_qp].addIa(_lambda * _volumetric_strain_increment.trace());
  _stress[_qp] += 2.0 * _G * _deviatoric_strain_increment - 2.0 * _G * delta_gamma * n;
}

void
Chaboche::updateIsotropicHardening(const Real eqvpstrain)
{
  // Isotropic hardening parameters
  Real Q = _Q->value(_t, _q_point[_qp]);
  Real b = _b->value(_t, _q_point[_qp]);
  Real sigma_0 = _sigma_0->value(_t, _q_point[_qp]);
	
  _isotropic_hardening[_qp] = sigma_0 + Q * (1.0 - std::exp(-b * eqvpstrain));
}

void
Chaboche::updateBackstress(const Real delta_gamma,
                           const RankTwoTensor n)
{
  _backstress1[_qp] += (2.0/3.0) * _C1->value(_t, _q_point[_qp]) * delta_gamma * n;
  _backstress1[_qp] -= _gamma1->value(_t, _q_point[_qp]) * delta_gamma * _backstress1[_qp];
  
  _backstress2[_qp] += (2.0/3.0) * _C2->value(_t, _q_point[_qp]) * delta_gamma * n;
  _backstress2[_qp] -= _gamma2->value(_t, _q_point[_qp]) * delta_gamma * _backstress2[_qp];
}
