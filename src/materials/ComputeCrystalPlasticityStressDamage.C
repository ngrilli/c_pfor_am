// Nicolò Grilli
// University of Bristol
// 17 Luglio 2022

#include "ComputeCrystalPlasticityStressDamage.h"

#include "CrystalPlasticityDislocationUpdateBase.h"
#include "libmesh/utility.h"
#include "Conversion.h"
#include "MooseException.h"

registerMooseObject("TensorMechanicsApp", ComputeCrystalPlasticityStressDamage);

InputParameters
ComputeCrystalPlasticityStressDamage::validParams()
{
  InputParameters params = ComputeFiniteStrainElasticStress::validParams();

  params.addClassDescription(
      "Crystal Plasticity base class: handles the Newton iteration over the stress residual and "
      "calculates the Jacobian based on constitutive laws from multiple material classes "
      "that are inherited from CrystalPlasticityDislocationUpdateBase."
	  "The difference between this class and ComputeMultipleCrystalPlasticityStress "
	  "is that the _models variable here is an array of CrystalPlasticityDislocationUpdateBase "
	  "instead of CrystalPlasticityStressUpdateBase. "
	  "This material model is coupled with phase field damage. ");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>(
      "E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<Real>("plastic_damage_prefactor",0.0,
                        "prefactor applied to the plastic work to determine the fraction "
						"of plastic energy that contributes to damage. ");
  params.addParam<std::string>(
      "base_name",
      "Optional parameter that allows the user to define multiple mechanics material systems on "
      "the same block, i.e. for multiple phases");

  params.addRequiredParam<std::vector<MaterialName>>(
      "crystal_plasticity_models",
      "The material objects to use to calculate crystal plasticity stress and strains.");
  params.addParam<std::vector<MaterialName>>("eigenstrain_names",
                                             "The material objects to calculate eigenstrains.");
  params.addParam<MooseEnum>("tan_mod_type",
                             MooseEnum("exact none", "none"),
                             "Type of tangent moduli for preconditioner: default elastic");
  params.addParam<Real>("rtol", 1e-6, "Constitutive stress residual relative tolerance");
  params.addParam<Real>("abs_tol", 1e-6, "Constitutive stress residual absolute tolerance");
  params.addParam<unsigned int>("maxiter", 100, "Maximum number of iterations for stress update");
  params.addParam<unsigned int>(
      "maxiter_state_variable", 100, "Maximum number of iterations for state variable update");
  params.addParam<unsigned int>(
      "maximum_substep_iteration", 1, "Maximum number of substep iteration");
  params.addParam<bool>("use_line_search", false, "Use line search in constitutive update");
  params.addParam<Real>("min_line_search_step_size", 0.01, "Minimum line search step size");
  params.addParam<Real>("line_search_tol", 0.5, "Line search bisection method tolerance");
  params.addParam<unsigned int>(
      "line_search_maxiter", 20, "Line search bisection method maximum number of iteration");
  params.addParam<MooseEnum>("line_search_method",
                             MooseEnum("CUT_HALF BISECTION", "CUT_HALF"),
                             "The method used in line search");
  params.addParam<bool>(
      "print_state_variable_convergence_error_messages",
      false,
      "Whether or not to print warning messages from the crystal plasticity specific convergence "
      "checks on the stress measure and general constitutive model quantinties.");
  params.addParam<UserObjectName>("read_initial_Fp",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element value "
                                  "of the initial plastic deformation gradient");
  return params;
}

ComputeCrystalPlasticityStressDamage::ComputeCrystalPlasticityStressDamage(
    const InputParameters & parameters)
  : ComputeFiniteStrainElasticStress(parameters),
    _num_models(getParam<std::vector<MaterialName>>("crystal_plasticity_models").size()),
    _num_eigenstrains(getParam<std::vector<MaterialName>>("eigenstrain_names").size()),
	_c(coupledValue("c")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _H(declareProperty<Real>("hist")), // History variable to avoid damage decrease 
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _E(declareProperty<Real>(getParam<MaterialPropertyName>("E_name"))),
    _dEdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"),
                                          getVar("c", 0)->name())),
    _d2Ed2c(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _dstress_dc(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())), 
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),		
	_D(getMaterialProperty<Real>("D_name")),
    _dDdc(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name())),
    _d2Dd2c(getMaterialPropertyDerivative<Real>(
        "D_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _fp_increment(declareProperty<RankTwoTensor>("Fp_increment")), // increment of Fp over _dt
	_plastic_work(declareProperty<Real>("plastic_work")), // scalar plastic work
	_plastic_work_old(getMaterialPropertyOld<Real>("plastic_work")), // and value at previous time step
	_elastic_deformation_grad(declareProperty<RankTwoTensor>("elastic_deformation_grad")), // Fe
    _plastic_damage_prefactor(getParam<Real>("plastic_damage_prefactor")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _rtol(getParam<Real>("rtol")),
    _abs_tol(getParam<Real>("abs_tol")),
    _maxiter(getParam<unsigned int>("maxiter")),
    _maxiterg(getParam<unsigned int>("maxiter_state_variable")),
    _tan_mod_type(getParam<MooseEnum>("tan_mod_type").getEnum<TangentModuliType>()),
    _max_substep_iter(getParam<unsigned int>("maximum_substep_iteration")),
    _use_line_search(getParam<bool>("use_line_search")),
    _min_line_search_step_size(getParam<Real>("min_line_search_step_size")),
    _line_search_tolerance(getParam<Real>("line_search_tol")),
    _line_search_max_iterations(getParam<unsigned int>("line_search_maxiter")),
    _line_search_method(getParam<MooseEnum>("line_search_method").getEnum<LineSearchMethod>()),
    _plastic_deformation_gradient(declareProperty<RankTwoTensor>("plastic_deformation_gradient")),
    _plastic_deformation_gradient_old(
        getMaterialPropertyOld<RankTwoTensor>("plastic_deformation_gradient")),
    _eigenstrain_deformation_gradient(
        _num_eigenstrains ? &declareProperty<RankTwoTensor>("eigenstrain_deformation_gradient")
                          : nullptr),
    _eigenstrain_deformation_gradient_old(
        _num_eigenstrains
            ? &getMaterialPropertyOld<RankTwoTensor>("eigenstrain_deformation_gradient")
            : nullptr),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _pk2(declareProperty<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _pk2_old(getMaterialPropertyOld<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _total_lagrangian_strain(
        declareProperty<RankTwoTensor>("total_lagrangian_strain")), // Lagrangian strain
    _updated_rotation(declareProperty<RankTwoTensor>("updated_rotation")),
    _crysrot(getMaterialProperty<RankTwoTensor>(
        "crysrot")), // defined in the elasticity tensor classes for crystal plasticity
    _print_convergence_message(getParam<bool>("print_state_variable_convergence_error_messages")),
    
	// UserObject to read the initial plastic deformation gradient from file						
    _read_initial_Fp(isParamValid("read_initial_Fp")
                               ? &getUserObject<ElementPropertyReadFile>("read_initial_Fp")
                               : nullptr)
{
  _convergence_failed = false;
}

void
ComputeCrystalPlasticityStressDamage::initQpStatefulProperties()
{
  // temporary variable to store the initial plastic deformation gradient
  // read from file and then assign it to _plastic_deformation_gradient
  RankTwoTensor initial_Fp;	
	
  // Initialize Fp
  if (_read_initial_Fp) { // Read initial plastic deformation gradient from file

    // The file will have one row for each element
    // each row will contain the components
    // Fp_{11} Fp_{12} Fp_{13} Fp_{21} Fp_{22} Fp_{23} Fp_{31} Fp_{32} Fp_{33} 
	
    for (unsigned int i = 0; i < 3; ++i) {
	  for (unsigned int j = 0; j < 3; ++j) {
        initial_Fp(i,j) = _read_initial_Fp->getData(_current_elem, 3*i+j);
	  }
	}
  
    _plastic_deformation_gradient[_qp] = initial_Fp;
  
  } else { // Initialize uniform plastic deformation gradient to identity
  
    _plastic_deformation_gradient[_qp].zero();
    _plastic_deformation_gradient[_qp].addIa(1.0);
  
  }

  if (_num_eigenstrains)
  {
    (*_eigenstrain_deformation_gradient)[_qp].zero();
    (*_eigenstrain_deformation_gradient)[_qp].addIa(1.0);
  }
  else
  {
    // set to identity if no eigenstrain is added
    _inverse_eigenstrain_deformation_grad.zero();
    _inverse_eigenstrain_deformation_grad.addIa(1.0);
  }

  _pk2[_qp].zero();

  _total_lagrangian_strain[_qp].zero();

  _updated_rotation[_qp].zero();
  _updated_rotation[_qp].addIa(1.0);

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->setQp(_qp);
    _models[i]->initQpStatefulProperties();
  }

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    _eigenstrains[i]->setQp(_qp);
    _eigenstrains[i]->initQpStatefulProperties();
  }
}

void
ComputeCrystalPlasticityStressDamage::initialSetup()
{
  // get crystal plasticity models
  std::vector<MaterialName> model_names =
      getParam<std::vector<MaterialName>>("crystal_plasticity_models");

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    CrystalPlasticityDislocationUpdateBase * model =
        dynamic_cast<CrystalPlasticityDislocationUpdateBase *>(&getMaterialByName(model_names[i]));

    if (model)
    {
      _models.push_back(model);
      // TODO: check to make sure that the material model is compatible with this class
    }
    else
      mooseError("Model " + model_names[i] +
                 " is not compatible with ComputeCrystalPlasticityStressDamage");
  }

  // get crystal plasticity eigenstrains
  std::vector<MaterialName> eigenstrain_names =
      getParam<std::vector<MaterialName>>("eigenstrain_names");

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    ComputeCrystalPlasticityEigenstrainBase * eigenstrain =
        dynamic_cast<ComputeCrystalPlasticityEigenstrainBase *>(
            &getMaterialByName(eigenstrain_names[i]));

    if (eigenstrain)
      _eigenstrains.push_back(eigenstrain);
    else
      mooseError("Eigenstrain" + eigenstrain_names[i] +
                 " is not compatible with ComputeCrystalPlasticityStressDamage");
  }
}

void
ComputeCrystalPlasticityStressDamage::computeQpStress()
{
  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->setQp(_qp);

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
    _eigenstrains[i]->setQp(_qp);

  updateStress(_stress[_qp], _Jacobian_mult[_qp]); // This is NOT the exact jacobian
}

void
ComputeCrystalPlasticityStressDamage::updateStress(RankTwoTensor & cauchy_stress,
                                                     RankFourTensor & jacobian_mult)
{
  // Does not support face/boundary material property calculation
  if (isBoundaryMaterial())
    return;

  // Initialize substepping variables
  unsigned int substep_iter = 1;
  unsigned int num_substep = 1;

  _temporary_deformation_gradient_old = _deformation_gradient_old[_qp];
  if (_temporary_deformation_gradient_old.det() == 0)
    _temporary_deformation_gradient_old.addIa(1.0);

  _delta_deformation_gradient = _deformation_gradient[_qp] - _temporary_deformation_gradient_old;

  // Loop through all models and calculate the schmid tensor for the current state of the crystal
  // lattice
  // Not sure if we should pass in the updated or the original rotation here
  // If not, then we should not need to compute the flow direction every iteration here
  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->calculateFlowDirection(_crysrot[_qp]);

  do
  {
    _convergence_failed = false;
    preSolveQp();

    _substep_dt = _dt / num_substep;
    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->setSubstepDt(_substep_dt);

    // calculate F^{eigen} only when we have eigenstrain
    if (_num_eigenstrains)
      calculateEigenstrainDeformationGrad();

    for (unsigned int istep = 0; istep < num_substep; ++istep)
    {
      _temporary_deformation_gradient =
          (static_cast<Real>(istep) + 1) / num_substep * _delta_deformation_gradient;
      _temporary_deformation_gradient += _temporary_deformation_gradient_old;

      solveQp();

      if (_convergence_failed)
      {
        if (_print_convergence_message)
          mooseWarning(
              "The crystal plasticity constitutive model has failed to converge. Increasing "
              "the number of substeps.");

        substep_iter++;
        num_substep *= 2;
        break;
      }
    }

    if (substep_iter > _max_substep_iter && _convergence_failed)
      mooseException("ComputeCrystalPlasticityStressDamage: Constitutive failure");
  } while (_convergence_failed);

  postSolveQp(cauchy_stress, jacobian_mult);
}

void
ComputeCrystalPlasticityStressDamage::preSolveQp()
{
  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->setInitialConstitutiveVariableValues();

  _pk2[_qp] = _pk2_old[_qp];
  _inverse_plastic_deformation_grad_old = _plastic_deformation_gradient_old[_qp].inverse();
}

void
ComputeCrystalPlasticityStressDamage::solveQp()
{
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->setSubstepConstitutiveVariableValues();
    _models[i]->calculateSlipResistance();
  }

  _inverse_plastic_deformation_grad = _inverse_plastic_deformation_grad_old;

  solveStateVariables();
  if (_convergence_failed)
    return; // pop back up and take a smaller substep

  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->updateSubstepConstitutiveVariableValues();

  // save off the old F^{p} inverse now that have converged on the stress and state variables
  _inverse_plastic_deformation_grad_old = _inverse_plastic_deformation_grad;
}

void
ComputeCrystalPlasticityStressDamage::postSolveQp(RankTwoTensor & cauchy_stress,
                                                    RankFourTensor & jacobian_mult)
{
  cauchy_stress = _elastic_deformation_gradient * _pk2[_qp] *
                  _elastic_deformation_gradient.transpose() / _elastic_deformation_gradient.det();

  calcTangentModuli(jacobian_mult);

  _total_lagrangian_strain[_qp] =
      _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] -
      RankTwoTensor::Identity();
  _total_lagrangian_strain[_qp] = _total_lagrangian_strain[_qp] * 0.5;

  // Calculate crystal rotation to track separately
  RankTwoTensor rot;
  _elastic_deformation_gradient.getRUDecompositionRotation(rot);
  _updated_rotation[_qp] = rot * _crysrot[_qp];
}

void
ComputeCrystalPlasticityStressDamage::solveStateVariables()
{
  unsigned int iteration;
  bool iter_flag = true;

  iteration = 0;
  // Check for slip system resistance update tolerance
  do
  {
    solveStress();
    if (_convergence_failed)
      return;

    _plastic_deformation_gradient[_qp] =
        _inverse_plastic_deformation_grad.inverse(); // the postSolveStress

    // Update slip system resistance and state variable after the stress has been finalized
    // We loop through all the models for each calculation
    // in order to make sure that when coupling appears, the state variables are updated based on
    // the same components
    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->cacheStateVariablesBeforeUpdate();

    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->calculateStateVariableEvolutionRateComponent();

    for (unsigned int i = 0; i < _num_models; ++i)
      if (!_models[i]->updateStateVariables())
        _convergence_failed = true;

    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->calculateSlipResistance();

    if (_convergence_failed)
      return;

    for (unsigned int i = 0; i < _num_models; ++i)
    {
      // iter_flag = false, stop iteration if all values are converged and good to go
      // iter_flag = true, continue iteration if any value is not converged
      if (!_models[i]->areConstitutiveStateVariablesConverged())
      {
        iter_flag = true;
        break;
      }
      else
        iter_flag = false; // iter_flag = false, stop iteration only when all models returns true
    }

    if (iter_flag)
    {
      if (_print_convergence_message)
        mooseWarning("ComputeCrystalPlasticityStressDamage: State variables (or the system "
                     "resistance) did not converge at element ",
                     _current_elem->id(),
                     " and qp ",
                     _qp,
                     "\n");
    }
    iteration++;
  } while (iter_flag && iteration < _maxiterg);

  if (iteration == _maxiterg)
  {
    if (_print_convergence_message)
      mooseWarning(
          "ComputeCrystalPlasticityStressDamage: Hardness Integration error. Reached the "
          "maximum number of iterations to solve for the state variables at element ",
          _current_elem->id(),
          " and qp ",
          _qp,
          "\n");

    _convergence_failed = true;
  }
}

void
ComputeCrystalPlasticityStressDamage::solveStress()
{
  unsigned int iteration = 0;
  RankTwoTensor dpk2;
  Real rnorm, rnorm0, rnorm_prev;

  // Calculate stress residual
  calculateResidualAndJacobian();
  if (_convergence_failed)
  {
    if (_print_convergence_message)
      mooseWarning("ComputeCrystalPlasticityStressDamage: the slip increment exceeds tolerance "
                   "at element ",
                   _current_elem->id(),
                   " and Gauss point ",
                   _qp);

    return;
  }

  rnorm = _residual_tensor.L2norm();
  rnorm0 = rnorm;

  // Check for stress residual tolerance; different from user object version which
  // compares the absolute tolerance of only the original rnorm value
  while (rnorm > _rtol * rnorm0 && rnorm > _abs_tol && iteration < _maxiter)
  {
    // Calculate stress increment
    dpk2 = -_jacobian.invSymm() * _residual_tensor;
    _pk2[_qp] = _pk2[_qp] + dpk2;

    calculateResidualAndJacobian();

    if (_convergence_failed)
    {
      if (_print_convergence_message)
        mooseWarning("ComputeCrystalPlasticityStressDamage: the slip increment exceeds tolerance "
                     "at element ",
                     _current_elem->id(),
                     " and Gauss point ",
                     _qp);

      return;
    }

    rnorm_prev = rnorm;
    rnorm = _residual_tensor.L2norm();

    if (_use_line_search && rnorm > rnorm_prev && !lineSearchUpdate(rnorm_prev, dpk2))
    {
      if (_print_convergence_message)
        mooseWarning("ComputeCrystalPlasticityStressDamage: Failed with line search");

      _convergence_failed = true;
      return;
    }

    if (_use_line_search)
      rnorm = _residual_tensor.L2norm();

    iteration++;
  }

  if (iteration >= _maxiter)
  {
    if (_print_convergence_message)
      mooseWarning("ComputeCrystalPlasticityStressDamage: Stress Integration error rmax = ",
                   rnorm,
                   " and the tolerance is ",
                   _rtol * rnorm0,
                   " when the rnorm0 value is ",
                   rnorm0,
                   " for element ",
                   _current_elem->id(),
                   " and qp ",
                   _qp);

    _convergence_failed = true;
  }
}

// Calculates stress residual equation and jacobian
void
ComputeCrystalPlasticityStressDamage::calculateResidualAndJacobian()
{
  calculateResidual();
  if (_convergence_failed)
    return;
  calculateJacobian();
}

void
ComputeCrystalPlasticityStressDamage::calculateResidual()
{
  RankTwoTensor ce, elastic_strain, ce_pk2, equivalent_slip_increment_per_model,
      equivalent_slip_increment, pk2_new;

  equivalent_slip_increment.zero();

  // calculate slip rate in order to compute F^{p-1}
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    equivalent_slip_increment_per_model.zero();

    // calculate shear stress with consideration of contribution from other physics
    _models[i]->calculateShearStress(
        _pk2[_qp], _inverse_eigenstrain_deformation_grad, _num_eigenstrains);

    _convergence_failed = !_models[i]->calculateSlipRate();

    if (_convergence_failed)
      return;

    _models[i]->calculateEquivalentSlipIncrement(equivalent_slip_increment_per_model);
    equivalent_slip_increment += equivalent_slip_increment_per_model;
  }

  RankTwoTensor residual_equivalent_slip_increment =
      RankTwoTensor::Identity() - equivalent_slip_increment;
  _inverse_plastic_deformation_grad =
      _inverse_plastic_deformation_grad_old * residual_equivalent_slip_increment;

  _elastic_deformation_gradient = _temporary_deformation_gradient *
                                  _inverse_eigenstrain_deformation_grad *
                                  _inverse_plastic_deformation_grad;

  ce = _elastic_deformation_gradient.transpose() * _elastic_deformation_gradient;

  elastic_strain = ce - RankTwoTensor::Identity();
  elastic_strain *= 0.5;

  pk2_new = _elasticity_tensor[_qp] * elastic_strain;
  _residual_tensor = _pk2[_qp] - pk2_new;
}

void
ComputeCrystalPlasticityStressDamage::calculateJacobian()
{
  // may not need to cache the dfpinvdpk2 here. need to double check
  RankFourTensor dfedfpinv, deedfe, dfpinvdpk2, dfpinvdpk2_per_model;

  RankTwoTensor ffeiginv = _temporary_deformation_gradient * _inverse_eigenstrain_deformation_grad;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        dfedfpinv(i, j, k, j) = ffeiginv(i, k);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->calculateTotalPlasticDeformationGradientDerivative(
        dfpinvdpk2_per_model,
        _inverse_plastic_deformation_grad_old,
        _inverse_eigenstrain_deformation_grad,
        _num_eigenstrains);
    dfpinvdpk2 += dfpinvdpk2_per_model;
  }

  _jacobian =
      RankFourTensor::IdentityFour() - (_elasticity_tensor[_qp] * deedfe * dfedfpinv * dfpinvdpk2);
}

void
ComputeCrystalPlasticityStressDamage::computeStrainVolumetric(Real & F_pos, Real & F_neg, 
                                                             RankTwoTensor & ee, RankTwoTensor & ce, 
														     RankTwoTensor & pk2_new)
{
  //Anisotropic elasticity
  Real Kb = 0.0; // Kb is the reference bulk modulus
  Real Je; // Je is relative elastic volume change
  Real Je23; // Je^{2/3}
  Real delta; // delta is the trace of volumetric part of elastic Green-Lagrange strain
  
  // Positive and negative volumetric free energies
  Real a_pos_vol, a_neg_vol;
  
  // Isochoric free energy
  Real a_pos_cpl;
  
  // Undamaged elastic energy tensor
  RankTwoTensor undamaged_elastic_energy;
  
  // Positive and negative part of the Piola-Kirchhoff stress
  RankTwoTensor pk2_pos;
  RankTwoTensor pk2_neg;
  
  // inverse of the right Cauchy–Green deformation tensor (elastic)
  RankTwoTensor invce;
  
  // Anisotropic elasticity (Luscher 2017)
  // Kb = K in (Luscher 2017)
  // Kb = (1/9) I : C : I
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      Kb += _elasticity_tensor[_qp](i, i, j, j);
  
  Kb = Kb / 9.0;

  Je = _elastic_deformation_gradient.det();
  Je23 = std::pow(Je,2.0/3.0);
  delta = 1.5 * (Je23 - 1.0);
  
  // Calculate volumetric free energy
  // Equations 15 and 16 in Grilli, Koslowski, 2019
  
  if (Je >= 1.0) { // expansion
	  
    a_pos_vol = 0.5 * Kb * delta * delta;
    a_neg_vol = 0.0;
	  
  } else { // compression
	  
    a_pos_vol = 0.0;
	a_neg_vol = 0.5 * Kb * delta * delta;
	
  }
  
  // Calculate isochoric free energy
  // Equation 17 in Grilli, Koslowski, 2019

  undamaged_elastic_energy = _elasticity_tensor[_qp] * ee;
  undamaged_elastic_energy = 0.5 * ee * undamaged_elastic_energy;

  a_pos_cpl = undamaged_elastic_energy.trace();
  a_pos_cpl -= 0.5 * Kb * delta * delta;

  // Calculate positive and negative parts of the Piola-Kirchhoff stress
  // Equation 18 in Grilli, Koslowski, 2019
  if (Je >= 1.0) { // expansion
  
    pk2_pos = _elasticity_tensor[_qp] * ee;
    pk2_neg = 0.0;
	
  } else { // compression
  
    invce = ce.inverse();
	
	pk2_neg = Je23 * Kb * delta * invce;
	
    pk2_pos = (-1.0) * pk2_neg;	
    pk2_pos += _elasticity_tensor[_qp] * ee;

  }

  // Positive part of the stress is degraded by _D[_qp]
  pk2_new = (_D[_qp] * pk2_pos) + pk2_neg;
  
  // Positive and negative parts of the free energy
  // Equations 13 and 14 in Grilli, Koslowski, 2019
  // additionally, plastic work is included for damage
  F_pos = a_pos_vol + a_pos_cpl + _plastic_damage_prefactor * _plastic_work[_qp];
  F_neg = a_neg_vol;
  
  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = pk2_pos * _dDdc[_qp];
  
  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history variable
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = pk2_pos * _dDdc[_qp];
  
}

void
ComputeCrystalPlasticityStressDamage::calcTangentModuli(RankFourTensor & jacobian_mult)
{
  switch (_tan_mod_type)
  {
    case TangentModuliType::EXACT:
      elastoPlasticTangentModuli(jacobian_mult);
      break;
    default:
      elasticTangentModuli(jacobian_mult);
  }
}

void
ComputeCrystalPlasticityStressDamage::elastoPlasticTangentModuli(RankFourTensor & jacobian_mult)
{
  RankFourTensor tan_mod;
  RankTwoTensor pk2fet, fepk2, feiginvfpinv;
  RankFourTensor deedfe, dsigdpk2dfe, dfedf;

  // Fill in the matrix stiffness material property
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  usingTensorIndices(i_, j_, k_, l_);
  dsigdpk2dfe = _elastic_deformation_gradient.times<i_, k_, j_, l_>(_elastic_deformation_gradient) *
                _elasticity_tensor[_qp] * deedfe;

  pk2fet = _pk2[_qp] * _elastic_deformation_gradient.transpose();
  fepk2 = _elastic_deformation_gradient * _pk2[_qp];

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
      {
        tan_mod(i, j, i, l) += pk2fet(l, j);
        tan_mod(i, j, j, l) += fepk2(i, l);
      }

  tan_mod += dsigdpk2dfe;

  const auto je = _elastic_deformation_gradient.det();
  if (je > 0.0)
    tan_mod /= je;

  feiginvfpinv = _inverse_eigenstrain_deformation_grad * _inverse_plastic_deformation_grad;
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
        dfedf(i, j, i, l) = feiginvfpinv(l, j);

  jacobian_mult = tan_mod * dfedf;
}

void
ComputeCrystalPlasticityStressDamage::elasticTangentModuli(RankFourTensor & jacobian_mult)
{
  // update jacobian_mult
  jacobian_mult = _elasticity_tensor[_qp];
}

bool
ComputeCrystalPlasticityStressDamage::lineSearchUpdate(const Real & rnorm_prev,
                                                         const RankTwoTensor & dpk2)
{
  if (_line_search_method == LineSearchMethod::CutHalf)
  {
    Real rnorm;
    Real step = 1.0;

    do
    {
      _pk2[_qp] = _pk2[_qp] - step * dpk2;
      step /= 2.0;
      _pk2[_qp] = _pk2[_qp] + step * dpk2;

      calculateResidual();
      rnorm = _residual_tensor.L2norm();
    } while (rnorm > rnorm_prev && step > _min_line_search_step_size);

    // has norm improved or is the step still above minumum search step size?
    return (rnorm <= rnorm_prev || step > _min_line_search_step_size);
  }
  else if (_line_search_method == LineSearchMethod::Bisection)
  {
    unsigned int count = 0;
    Real step_a = 0.0;
    Real step_b = 1.0;
    Real step = 1.0;
    Real s_m = 1000.0;
    Real rnorm = 1000.0;

    calculateResidual();
    auto s_b = _residual_tensor.doubleContraction(dpk2);
    const auto rnorm1 = _residual_tensor.L2norm();
    _pk2[_qp] = _pk2[_qp] - dpk2;
    calculateResidual();
    auto s_a = _residual_tensor.doubleContraction(dpk2);
    const auto rnorm0 = _residual_tensor.L2norm();
    _pk2[_qp] = _pk2[_qp] + dpk2;

    if ((rnorm1 / rnorm0) < _line_search_tolerance || s_a * s_b > 0)
    {
      calculateResidual();
      return true;
    }

    while ((rnorm / rnorm0) > _line_search_tolerance && count < _line_search_max_iterations)
    {
      _pk2[_qp] = _pk2[_qp] - step * dpk2;
      step = 0.5 * (step_b + step_a);
      _pk2[_qp] = _pk2[_qp] + step * dpk2;
      calculateResidual();
      s_m = _residual_tensor.doubleContraction(dpk2);
      rnorm = _residual_tensor.L2norm();

      if (s_m * s_a < 0.0)
      {
        step_b = step;
        s_b = s_m;
      }
      if (s_m * s_b < 0.0)
      {
        step_a = step;
        s_a = s_m;
      }
      count++;
    }

    // below tolerance and max iterations?
    return ((rnorm / rnorm0) < _line_search_tolerance && count < _line_search_max_iterations);
  }
  else
    mooseError("Line search method is not provided.");
}

void
ComputeCrystalPlasticityStressDamage::calculateEigenstrainDeformationGrad()
{
  _inverse_eigenstrain_deformation_grad.zero();
  _inverse_eigenstrain_deformation_grad.addIa(1.0);

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    _eigenstrains[i]->setSubstepDt(_substep_dt);
    _eigenstrains[i]->computeQpProperties();
    _inverse_eigenstrain_deformation_grad *= _eigenstrains[i]->getDeformationGradientInverse();
  }
  (*_eigenstrain_deformation_gradient)[_qp] = _inverse_eigenstrain_deformation_grad.inverse();
}
