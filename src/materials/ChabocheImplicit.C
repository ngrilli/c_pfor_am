// Nicolò Grilli
// Università di Bristol
// 28 Maggio 2026

#include "ChabocheImplicit.h"

registerMooseObject("c_pfor_amApp", ChabocheImplicit);

InputParameters
ChabocheImplicit::validParams()
{
  InputParameters params = ComputeStressBase::validParams();

  params.addClassDescription("A Chaboche model with return mapping, implicit version");
  params.addParam<MaterialPropertyName>("sigma_0_name","sigma_0","Constant stress of the isotropic hardening");
  params.addParam<MaterialPropertyName>("Q_name","Q","Maximum isotropic hardening");
  params.addParam<MaterialPropertyName>("b_name","b","Hardening rate");
  params.addParam<MaterialPropertyName>("E_name","E","Young's modulus");
  params.addParam<MaterialPropertyName>("nu_name","nu","Poisson's ratio");
  params.addParam<MaterialPropertyName>("C1_name","C1","C constant for the first backstress");
  params.addParam<MaterialPropertyName>("gamma1_name","gamma1","gamma constant for the first backstress");
  params.addParam<MaterialPropertyName>("C2_name","C2","C constant for the second backstress");
  params.addParam<MaterialPropertyName>("gamma2_name","gamma2","gamma constant for the second backstress");
  params.addParam<Real>("tolerance",1e-6,"Yield function tolerance");
  params.addParam<int>("max_iterations",1000,"Number of return mapping iterations before unconverged");
  params.addParam<MooseEnum>("tangent_moduli_type",
                             MooseEnum("elastic elasto_plastic", "elastic"),
                             "Type of tangent moduli for preconditioner: default elastic");
  return params;
}

ChabocheImplicit::ChabocheImplicit(const InputParameters & parameters)
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
    _sigma_0(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("sigma_0_name"))),
    _Q(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("Q_name"))),
    _b(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("b_name"))),
    _E(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("E_name"))),
    _nu(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("nu_name"))),
    _C1(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("C1_name"))),
    _gamma1(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("gamma1_name"))),
    _C2(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("C2_name"))),
    _gamma2(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("gamma2_name"))),
    _tolerance(getParam<Real>("tolerance")),
    _max_iterations(getParam<int>("max_iterations")),
    _tan_mod_type(getParam<MooseEnum>("tangent_moduli_type").getEnum<TangentModuliType>())
{
}

void
ChabocheImplicit::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  _plastic_strain[_qp].zero();
  _eqv_plastic_strain[_qp] = 0.0;
  _backstress1[_qp].zero();
  _backstress2[_qp].zero();
  _isotropic_hardening[_qp] = _sigma_0[_qp];
  initMapVoigt();
}

void
ChabocheImplicit::computeQpStress()
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
    _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
  }

  // Rotate the stress tensor to the current configuration
  _stress[_qp] = _rotation_increment[_qp] * _stress[_qp] * _rotation_increment[_qp].transpose();

  // Rotate plastic strain tensor to the current configuration
  _plastic_strain[_qp] =
      _rotation_increment[_qp] * _plastic_strain[_qp] * _rotation_increment[_qp].transpose();

  // Calculate the elastic strain
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _plastic_strain[_qp];
}

// Compute shear and bulk modulus
void
ChabocheImplicit::computeElasticConstants()
{
  _G = _E[_qp] / 2.0 / (1.0 + _nu[_qp]);
  _K = _E[_qp] / 3.0 / (1.0 - 2.0 * _nu[_qp]);
  _lambda = _K - (2.0 / 3.0) * _G;
}

void
ChabocheImplicit::decomposeStress(const RankTwoTensor & stress_old)
{
  _deviatoric_stress_old = stress_old.deviatoric();
  _volumetric_stress_old.zero();
  _volumetric_stress_old.addIa(stress_old.trace() / 3.0);
}

void
ChabocheImplicit::decomposeStrainIncrement(const RankTwoTensor & strain_increment)
{
  _deviatoric_strain_increment = strain_increment.deviatoric();
  _volumetric_strain_increment.zero();
  _volumetric_strain_increment.addIa(strain_increment.trace() / 3.0);
}

// Trial stress assuming the entire strain increment is elastic
// Trial stress includes only the deviatoric component
RankTwoTensor
ChabocheImplicit::computeTrialStress()
{
  return _deviatoric_stress_old + 2.0 * _G * _deviatoric_strain_increment;
}

// Calculate Mises equivalent stress
// this calculation intrinsically removes volumetric components
Real
ChabocheImplicit::getMisesEquivalent(const RankTwoTensor & stress)
{
  return std::sqrt(3 * stress.secondInvariant());
}

// Yield function, plasticity if positive
Real
ChabocheImplicit::yieldFunction(const RankTwoTensor & effective_deviatoric_stress,
                                const Real yield_stress)
{
  return getMisesEquivalent(effective_deviatoric_stress) - yield_stress;
}

void
ChabocheImplicit::returnMap(const Real eqvpstrain_old,
                            const RankTwoTensor & plastic_strain_old,
                            Real & eqvpstrain,
                            RankTwoTensor & plastic_strain)
{
  // scalar plastic multiplier that is solved for
  Real delta_gamma = 0.0;
  
  // Residual and Jacobian of the delta_gamma variable
  Real residual_delta_gamma;
  Real jacobian_delta_gamma;

  // Residual and Jacobian of the backstresses
  RankTwoTensor residual_backstress1;
  RankFourTensor jacobian_backstress1;
  RankTwoTensor residual_backstress2;
  RankFourTensor jacobian_backstress2;

  // off diagonal Jacobian terms
  RankTwoTensor dresidual_delta_gamma_dX; // derivative of the delta_gamma residual with respect to the backstresses
  RankTwoTensor dresidual_backstress1_ddelta_gamma; // derivative of the backstress1 residual with respect to delta_gamma
  RankTwoTensor dresidual_backstress2_ddelta_gamma; // derivative of the backstress2 residual with respect to delta_gamma
  RankFourTensor dresidual_backstress1_dbackstress2; // derivative of the backstress1 residual with respect to backstress2
  RankFourTensor dresidual_backstress2_dbackstress1; // derivative of the backstress2 residual with respect to backstress1

  // Direction of plastic flow, in tensor form, and its derivative with respect to backstresses
  RankTwoTensor n;
  RankFourTensor dn_dbackstress1;
  RankFourTensor dn_dbackstress2;
  
  // Scalar plastic strain increment
  Real delta_eps_p;
  
  // Mises equivalent of the effective deviatoric stress
  Real eqv_effective_deviatoric_stress;
  
  // Isotropic hardening parameters
  Real Q = _Q[_qp];
  Real b = _b[_qp];

  // I tensor product I = delta_{ik} delta_{jl}
  RankFourTensor I4(RankFourTensor::initIdentityFour);

  // symmetric part of I4 = 0.5*(delta_{ik} delta_{jl} + delta_{il} delta_{jk})
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  // full Jacobian for return mapping
  // 1 (delta_gamma) + 6 (symmetric backstress1) + 6 (symmetric backstress2) = 13
  std::vector<std::vector<Real>> J(13);
  for (auto & row : J)
    row.resize(13);
  // access components as J[row][col]

  // Initialize temporary backstress variables for the return mapping iterations
  _backstress1_iter = _backstress1[_qp];
  _backstress2_iter = _backstress2[_qp];

  // Initialize direction of plastic flow
  eqv_effective_deviatoric_stress = getMisesEquivalent(_effective_deviatoric_stress);
  n = (3.0 / 2.0) * _effective_deviatoric_stress / eqv_effective_deviatoric_stress;

  // Iteration counter
  unsigned int i = 0;
  
  // Staggered approach: first update delta_gamma at fixed backstresses, then update backstresses with the updated delta_gamma, repeat until convergence
  do {
    // Update plastic strain increment
    delta_eps_p = std::abs(delta_gamma);    
    
    // Update isotropic hardening
    eqvpstrain = eqvpstrain_old + delta_eps_p;
    _isotropic_hardening[_qp] = updateIsotropicHardening(eqvpstrain);

    // Update effective deviatoric stress, note that 
    // _trial_stress - _backstress1_iter - _backstress2_iter - 2.0 * _G * delta_gamma * n
    // is the same direction as
    // _trial_stress - _backstress1_iter - _backstress2_iter
    // because of Drucker principle of normality
    _effective_deviatoric_stress = _trial_stress - _backstress1_iter - _backstress2_iter;
    eqv_effective_deviatoric_stress = getMisesEquivalent(_effective_deviatoric_stress);
    n = (3.0 / 2.0) * _effective_deviatoric_stress / eqv_effective_deviatoric_stress;
    dresidual_delta_gamma_dX = -n; // same for both backstresses
    
    // Calculate delta_gamma residual
    residual_delta_gamma = eqv_effective_deviatoric_stress - _isotropic_hardening[_qp];

    // Calculate delta_gamma Jacobian
    jacobian_delta_gamma = -3 * _G - Q * b * std::exp(-b * eqvpstrain);

    // Update plastic multiplier
    delta_gamma -= residual_delta_gamma / jacobian_delta_gamma;

    // Calculate backstress residuals
    residual_backstress1 = _backstress1_iter - _backstress1_old[_qp] + _gamma1[_qp] * delta_gamma * _backstress1_iter - (2.0 / 3.0) * _C1[_qp] * delta_gamma * n;
    residual_backstress2 = _backstress2_iter - _backstress2_old[_qp] + _gamma2[_qp] * delta_gamma * _backstress2_iter - (2.0 / 3.0) * _C2[_qp] * delta_gamma * n;
    
    // Derivatives of the direction of plastic flow with respect to backstresses
    dn_dbackstress1 = ((-3.0 / 2.0) * I4 + n.outerProduct(n)) / eqv_effective_deviatoric_stress;
    dn_dbackstress2 = ((-3.0 / 2.0) * I4 + n.outerProduct(n)) / eqv_effective_deviatoric_stress;

    // Calculate backstress Jacobians
    jacobian_backstress1 = (1.0 + _gamma1[_qp] * delta_gamma) * I4 - (2.0 / 3.0) * _C1[_qp] * delta_gamma * dn_dbackstress1;
    jacobian_backstress2 = (1.0 + _gamma2[_qp] * delta_gamma) * I4 - (2.0 / 3.0) * _C2[_qp] * delta_gamma * dn_dbackstress2;

    // off diagonal Jacobian of the backstress residuals with respect to delta_gamma
    dresidual_backstress1_ddelta_gamma = _gamma1[_qp] * _backstress1_iter - (2.0 / 3.0) * _C1[_qp] * n;
    dresidual_backstress2_ddelta_gamma = _gamma2[_qp] * _backstress2_iter - (2.0 / 3.0) * _C2[_qp] * n;

    // off diagonal Jacobians of the backstress residuals with respect to backstresses
    dresidual_backstress1_dbackstress2 = - (2.0 / 3.0) * _C1[_qp] * delta_gamma * dn_dbackstress2;
    dresidual_backstress2_dbackstress1 = - (2.0 / 3.0) * _C2[_qp] * delta_gamma * dn_dbackstress1;

    // Update backstresses with a Newton step
    _backstress1_iter -= jacobian_backstress1.inverse() * residual_backstress1;
    _backstress2_iter -= jacobian_backstress2.inverse() * residual_backstress2;

    if (i > _max_iterations) // unconverged
      mooseError("Constitutive failure");
    
    // update iteration counter
    i++;
  }
  while (residual_delta_gamma > _tolerance || residual_backstress1.L2norm() > _tolerance || residual_backstress2.L2norm() > _tolerance);
  
  // Update plastic strain, isotropic hardening and back stress
  delta_eps_p = std::abs(delta_gamma);
  eqvpstrain = eqvpstrain_old + delta_eps_p;
  _isotropic_hardening[_qp] = updateIsotropicHardening(eqvpstrain);

  _effective_deviatoric_stress = _trial_stress - _backstress1_iter - _backstress2_iter;
  eqv_effective_deviatoric_stress = getMisesEquivalent(_effective_deviatoric_stress);
  n = (3.0 / 2.0) * _effective_deviatoric_stress / eqv_effective_deviatoric_stress;  
  plastic_strain = plastic_strain_old + delta_gamma * n;

  _backstress1[_qp] = _backstress1_iter;
  _backstress2[_qp] = _backstress2_iter;

  // Update total stress: recombine deviatoric and hydrostatic parts
  _stress[_qp] = _stress_old[_qp] + 2.0 * _G * _volumetric_strain_increment;
  _stress[_qp].addIa(_lambda * _volumetric_strain_increment.trace());
  _stress[_qp] += 2.0 * _G * _deviatoric_strain_increment - 2.0 * _G * delta_gamma * n;

  switch (_tan_mod_type)
  {
    case TangentModuliType::ELASTO_PLASTIC:
      elastoPlasticTangentModuli(eqvpstrain, n);
      break;
    default:
      _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
  }
}

Real
ChabocheImplicit::updateIsotropicHardening(const Real eqvpstrain)
{
  // Isotropic hardening parameters
  Real Q = _Q[_qp];
  Real b = _b[_qp];
  Real sigma_0 = _sigma_0[_qp];
	
  return sigma_0 + Q * (1.0 - std::exp(-b * eqvpstrain));
}

RankTwoTensor
ChabocheImplicit::updateBackstress(const Real delta_gamma,
                                   const RankTwoTensor n,
                                   const unsigned int backstress_index)
{
  if (backstress_index == 0) {
    return (_backstress1_old[_qp] + (2.0 / 3.0) * _C1[_qp] * delta_gamma * n) /
           (1.0 + _gamma1[_qp] * delta_gamma);
  } else if (backstress_index == 1) {
    return (_backstress2_old[_qp] + (2.0 / 3.0) * _C2[_qp] * delta_gamma * n) /
           (1.0 + _gamma2[_qp] * delta_gamma);
  } else {
    mooseError("Invalid backstress index");
  }
}

Real
ChabocheImplicit::numericalJacobian(const Real delta_gamma,
                                    const Real eqvpstrain_old,
                                    const RankTwoTensor n,
                                    const Real residual)
{
  // Perturb delta_gamma
  Real h = 1e-8;
  Real delta_gamma_perturbed = delta_gamma + h;

  Real delta_eps_p_perturbed = std::abs(delta_gamma_perturbed);
  Real eqvpstrain_perturbed = eqvpstrain_old + delta_eps_p_perturbed;

  // Update isotropic hardening with perturbed delta_gamma
  Real isotropic_hardening_perturbed = updateIsotropicHardening(eqvpstrain_perturbed);

  // Update backstress with perturbed delta_gamma
  RankTwoTensor backstress1_perturbed = updateBackstress(delta_gamma_perturbed, n, 0);
  RankTwoTensor backstress2_perturbed = updateBackstress(delta_gamma_perturbed, n, 1);

  // Update effective deviatoric stress with perturbed delta_gamma
  RankTwoTensor effective_deviatoric_stress_perturbed = _trial_stress - backstress1_perturbed - backstress2_perturbed - 2.0 * _G * delta_gamma_perturbed * n;
  Real eqv_effective_deviatoric_stress_perturbed = getMisesEquivalent(effective_deviatoric_stress_perturbed);

  // Calculate residual with perturbed delta_gamma
  Real residual_perturbed = eqv_effective_deviatoric_stress_perturbed - isotropic_hardening_perturbed;

  // Calculate numerical Jacobian
  return (residual_perturbed - residual) / h;
}

void
ChabocheImplicit::elastoPlasticTangentModuli(const Real eqvpstrain,
                                             const RankTwoTensor n)
{
  RankFourTensor tan_mod;
  RankFourTensor I4(RankFourTensor::initIdentityFour); // I tensor product I
  //RankFourTensor I4dev(RankFourTensor::initIdentityDeviatoric); // deviatoric 4th order identity
  RankFourTensor I4dev(RankFourTensor::initIdentitySymmetricFour);
  RankFourTensor n_outer_n; // tensor product of n with itself

  // Hardening rate
  Real H_iso = _Q[_qp] * _b[_qp] * std::exp(-_b[_qp] * eqvpstrain);
  Real H_kin = (2.0 / 3.0) * (_C1[_qp] + _C2[_qp]);
  Real H = H_iso + H_kin;

  n_outer_n = n.outerProduct(n);
  
  _Jacobian_mult[_qp] = _K * I4 + 2.0 * _G * (1.0 - (2.0 * _G) / (2.0 * _G + H)) * (I4dev - (3.0/2.0) * n_outer_n);
}

// Convert symmetric 3x3x3x3 matrix into 6x6 matrix using Mandel notation
std::vector<std::vector<Real>>
ChabocheImplicit::convertSym3333ToMandel66(const RankFourTensor tensor)
{  
  std::vector<std::vector<Real>> matrix(6);
  for (auto & row : matrix)
    row.resize(6);

    // TO DO

  return matrix;
}

// Map from symmetric 3x3 tensor indices to Mandel-Voigt 6x1 vector indices
void
ChabocheImplicit::initMapVoigt()
{
  _map_Voigt.resize(2);
  for (auto & row : _map_Voigt)
    row.resize(6);

  _map_Voigt[0][0] = 0; // xx
  _map_Voigt[1][0] = 0;

  _map_Voigt[0][1] = 1; // yy
  _map_Voigt[1][1] = 1;

  _map_Voigt[0][2] = 2; // zz
  _map_Voigt[1][2] = 2;

  _map_Voigt[0][3] = 1; // yz
  _map_Voigt[1][3] = 2;

  _map_Voigt[0][4] = 0; // xz
  _map_Voigt[1][4] = 2;

  _map_Voigt[0][5] = 0; // xy
  _map_Voigt[1][5] = 1;

  _weight_Mandel.resize(6);
  _weight_Mandel[0] = 1.0; // xx
  _weight_Mandel[1] = 1.0; // yy
  _weight_Mandel[2] = 1.0; // zz
  _weight_Mandel[3] = std::sqrt(2.0); // yz
  _weight_Mandel[4] = std::sqrt(2.0); // xz
  _weight_Mandel[5] = std::sqrt(2.0); // xy
}
