// Nicolò Grilli
// University of Bristol
// 15 Luglio 2023

#include "CrystalPlasticityUndamagedStress.h"

#include "CrystalPlasticityDislocationUpdateBase.h"
#include "ComputeCrystalPlasticityEigenstrainBase.h"
#include "libmesh/utility.h"
#include "Conversion.h"
#include "MooseException.h"

registerMooseObject("TensorMechanicsApp", CrystalPlasticityUndamagedStress);

InputParameters
CrystalPlasticityUndamagedStress::validParams()
{
  InputParameters params = ComputeMultipleCrystalPlasticityStress::validParams();
  params.addClassDescription(
      "This is very similar to ComputeCrystalPlasticityStressDamage "
      "but the plastic strain rate on the slip systems is calculated with the undamaged stress "
      "while the damaged stress and jacobian is passed to stress equilibrium calculations. ");
  params.addCoupledVar("c", 0.0, "Order parameter for damage");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>(
      "E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<Real>("plastic_damage_prefactor", 0.0,
                        "prefactor applied to the plastic work to determine the fraction "
  						"of plastic energy that contributes to damage. ");
  params.addParam<UserObjectName>("read_initial_Fp",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element value "
                                  "of the initial plastic deformation gradient");
  params.addCoupledVar("temperature",303.0,"Temperature");
  params.addParam<Real>("reference_temperature",303.0,"reference temperature for thermal expansion");
  params.addParam<Real>("thermal_expansion",0.0,"Linear thermal expansion coefficient");
  params.addParam<Real>("dCTE_dT",0.0,"First derivative of the thermal expansion coefficient with respect to temperature");
  params.addParam<bool>("suppress_constitutive_failure", false, "Use old values of Fp if NR algorithm fails.");
  params.addParam<bool>("use_snes_vi_solver",
                        false,
                        "Use PETSc's SNES variational inequalities solver to enforce damage "
                        "irreversibility condition and restrict damage value <= 1.");
  return params;
}

CrystalPlasticityUndamagedStress::CrystalPlasticityUndamagedStress(
    const InputParameters & parameters)
  : ComputeMultipleCrystalPlasticityStress(parameters),
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
    _plastic_damage_prefactor(getParam<Real>("plastic_damage_prefactor")),

	// UserObject to read the initial plastic deformation gradient from file						
    _read_initial_Fp(isParamValid("read_initial_Fp")
                               ? &getUserObject<ElementPropertyReadFile>("read_initial_Fp")
                               : nullptr),
                               
    // Thermal expansion variables
    _temperature(coupledValue("temperature")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _dCTE_dT(getParam<Real>("dCTE_dT")),
    
    _suppress_constitutive_failure(getParam<bool>("suppress_constitutive_failure")),
    _use_snes_vi_solver(getParam<bool>("use_snes_vi_solver")),
    
    // _pk2 is the undamaged stress used for the crystal plasticity NR algorithm
    _pk2_damaged(declareProperty<RankTwoTensor>("pk2_damaged"))
{
  _convergence_failed = false;
}

void
CrystalPlasticityUndamagedStress::initQpStatefulProperties()
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
  _pk2_damaged[_qp].zero();

  _total_lagrangian_strain[_qp].zero();

  _updated_rotation[_qp].zero();
  _updated_rotation[_qp].addIa(1.0);

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _dislocation_models[i]->setQp(_qp);
    _dislocation_models[i]->initQpStatefulProperties();
  }

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    _eigenstrains[i]->setQp(_qp);
    _eigenstrains[i]->initQpStatefulProperties();
  }

  // Initialize plastic work
  _plastic_work[_qp] = 0.0;
  
  // Initialize increment of plastic deformation gradient
  _fp_increment[_qp].zero();
  
  // Initialize history variable
  _H[_qp] = 0.0;
}

void
CrystalPlasticityUndamagedStress::initialSetup()
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
      _dislocation_models.push_back(model);
      // TODO: check to make sure that the material model is compatible with this class
    }
    else
      mooseError("Model " + model_names[i] +
                 " is not compatible with CrystalPlasticityUndamagedStress");
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
                 " is not compatible with CrystalPlasticityUndamagedStress");
  }
}

void
CrystalPlasticityUndamagedStress::computeQpStress()
{
  for (unsigned int i = 0; i < _num_models; ++i)
    _dislocation_models[i]->setQp(_qp);

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
    _eigenstrains[i]->setQp(_qp);

  updateStress(_stress[_qp], _Jacobian_mult[_qp]);
}

void
CrystalPlasticityUndamagedStress::updateStress(RankTwoTensor & cauchy_stress,
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
    _dislocation_models[i]->calculateFlowDirection(_crysrot[_qp]);

  do
  {
    _convergence_failed = false;
    preSolveQp();

    _substep_dt = _dt / num_substep;
    for (unsigned int i = 0; i < _num_models; ++i)
      _dislocation_models[i]->setSubstepDt(_substep_dt);

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

    if (substep_iter > _max_substep_iter && _convergence_failed) {
		
      if (_suppress_constitutive_failure) {
		  
        _plastic_deformation_gradient[_qp] = _plastic_deformation_gradient_old[_qp];
        _convergence_failed = false;
		  
	  } else {
		  
        mooseException("CrystalPlasticityUndamagedStress: Constitutive failure");
		  
	  }
	}
      
  } while (_convergence_failed);

  postSolveQp(cauchy_stress, jacobian_mult);
}

void
CrystalPlasticityUndamagedStress::preSolveQp()
{
  for (unsigned int i = 0; i < _num_models; ++i)
    _dislocation_models[i]->setInitialConstitutiveVariableValues();

  _pk2[_qp] = _pk2_old[_qp];
  _inverse_plastic_deformation_grad_old = _plastic_deformation_gradient_old[_qp].inverse();
}

void
CrystalPlasticityUndamagedStress::solveQp()
{
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _dislocation_models[i]->setSubstepConstitutiveVariableValues();
    _dislocation_models[i]->calculateSlipResistance();
  }

  _inverse_plastic_deformation_grad = _inverse_plastic_deformation_grad_old;

  solveStateVariables();
  if (_convergence_failed)
    return; // pop back up and take a smaller substep

  for (unsigned int i = 0; i < _num_models; ++i)
    _dislocation_models[i]->updateSubstepConstitutiveVariableValues();

  // save off the old F^{p} inverse now that have converged on the stress and state variables
  _inverse_plastic_deformation_grad_old = _inverse_plastic_deformation_grad;
  
}

// cauchy_stress is calculated and passed to the stress equilibrium
// therefore it must be the damaged stress
// same for jacobian_mult
void
CrystalPlasticityUndamagedStress::postSolveQp(RankTwoTensor & cauchy_stress,
                                                    RankFourTensor & jacobian_mult)
{
  cauchy_stress = _elastic_deformation_gradient * _pk2_damaged[_qp] *
                  _elastic_deformation_gradient.transpose() / _elastic_deformation_gradient.det();

  calcTangentModuli(jacobian_mult);
  
  // Calculate increment of Fp over _dt
  _fp_increment[_qp] = _plastic_deformation_gradient[_qp] - _plastic_deformation_gradient_old[_qp];
  
  // update plastic work
  updatePlasticWork();

  _total_lagrangian_strain[_qp] =
      _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] -
      RankTwoTensor::Identity();
  _total_lagrangian_strain[_qp] = _total_lagrangian_strain[_qp] * 0.5;

  // Calculate crystal rotation to track separately
  RankTwoTensor rot;
  _elastic_deformation_gradient.getRUDecompositionRotation(rot);
  _updated_rotation[_qp] = rot * _crysrot[_qp];
}

// Plastic work updated according to
// equation (17) in:
// Elastic plastic deformation at finite strains
// E. H. Lee 1968,
// Stanford University technical report AD678483
void
CrystalPlasticityUndamagedStress::updatePlasticWork()
{
  // Temporary variable to store the tensor plastic work
  // the trace of this tensor is the scalar plastic work
  RankTwoTensor plastic_work_rate;
  
  Real Je; // Je is relative elastic volume change
  
  Je = _elastic_deformation_gradient.det();
  
  // Plastic work calculation is also based on damaged stress
  plastic_work_rate = Je * _stress[_qp] * _elastic_deformation_gradient 
  * _fp_increment[_qp] * _inverse_plastic_deformation_grad * _elastic_deformation_gradient.inverse();
	
  _plastic_work[_qp] = _plastic_work_old[_qp] + std::abs(plastic_work_rate.trace());
}

void
CrystalPlasticityUndamagedStress::solveStateVariables()
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
      _dislocation_models[i]->cacheStateVariablesBeforeUpdate();

    for (unsigned int i = 0; i < _num_models; ++i)
      _dislocation_models[i]->calculateStateVariableEvolutionRateComponent();

    for (unsigned int i = 0; i < _num_models; ++i)
      if (!_dislocation_models[i]->updateStateVariables())
        _convergence_failed = true;

    for (unsigned int i = 0; i < _num_models; ++i)
      _dislocation_models[i]->calculateSlipResistance();

    if (_convergence_failed)
      return;

    for (unsigned int i = 0; i < _num_models; ++i)
    {
      // iter_flag = false, stop iteration if all values are converged and good to go
      // iter_flag = true, continue iteration if any value is not converged
      if (!_dislocation_models[i]->areConstitutiveStateVariablesConverged())
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
        mooseWarning("CrystalPlasticityUndamagedStress: State variables (or the system "
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
          "CrystalPlasticityUndamagedStress: Hardness Integration error. Reached the "
          "maximum number of iterations to solve for the state variables at element ",
          _current_elem->id(),
          " and qp ",
          _qp,
          "\n");

    _convergence_failed = true;
  }
}

void
CrystalPlasticityUndamagedStress::solveStress()
{
  unsigned int iteration = 0;
  RankTwoTensor dpk2;
  Real rnorm, rnorm0, rnorm_prev;

  // Calculate stress residual
  calculateResidualAndJacobian();

  if (_convergence_failed)
  {
    if (_print_convergence_message)
      mooseWarning("CrystalPlasticityUndamagedStress: the slip increment exceeds tolerance "
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
    // NR algorithm is based on the undamaged stress
    dpk2 = -_jacobian.invSymm() * _residual_tensor;
    _pk2[_qp] = _pk2[_qp] + dpk2;

    calculateResidualAndJacobian();

    if (_convergence_failed)
    {
      if (_print_convergence_message)
        mooseWarning("CrystalPlasticityUndamagedStress: the slip increment exceeds tolerance "
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
        mooseWarning("CrystalPlasticityUndamagedStress: Failed with line search");

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
      mooseWarning("CrystalPlasticityUndamagedStress: Stress Integration error rmax = ",
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
CrystalPlasticityUndamagedStress::calculateResidualAndJacobian()
{
  calculateResidual();
  if (_convergence_failed)
    return;
  calculateJacobian();
}

void
CrystalPlasticityUndamagedStress::calculateResidual()
{
  RankTwoTensor ce, elastic_strain, equivalent_slip_increment_per_model,
      equivalent_slip_increment, pk2_new;
      
  Real F_pos, F_neg; // tensile and compressive part of the elastic strain energy
  
  // Thermal expansion variables
  Real temperature = _temperature[_qp];
  Real reference_temperature = _reference_temperature;
  Real thermal_expansion = _thermal_expansion;
  Real dCTE_dT =_dCTE_dT;

  equivalent_slip_increment.zero();

  // calculate slip rate in order to compute F^{p-1}
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    equivalent_slip_increment_per_model.zero();

    // calculate shear stress with consideration of contribution from other physics
    // resolved shear stress is calculated with the undamaged stress
    // because it is used to calculate the plastic strain rate
    _dislocation_models[i]->calculateShearStress(
        _pk2[_qp], _inverse_eigenstrain_deformation_grad, _num_eigenstrains);

    _convergence_failed = !_dislocation_models[i]->calculateSlipRate();

    if (_convergence_failed)
      return;

    _dislocation_models[i]->calculateEquivalentSlipIncrement(equivalent_slip_increment_per_model);
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
  
  // Calculate volumetric thermal expansion and thermal eigenstrain
  _volumetric_thermal_expansion = 1.5 * (
    std::exp((2.0/3.0)*(
        0.5 * dCTE_dT * (temperature-reference_temperature)*(temperature-reference_temperature)
        + thermal_expansion * (temperature - reference_temperature)
      )
    ) - 1.0);

  _thermal_eigenstrain = (1.0/3.0) * _volumetric_thermal_expansion * RankTwoTensor::Identity();
  
  // Decompose ee into volumetric and non-volumetric
  // and calculate elastic energy and stress
  // pk2_new is the undamaged stress
  computeStrainVolumetric(F_pos, F_neg, elastic_strain, ce, pk2_new);
  
  // calculate history variable and
  // assign elastic free energy to _E
  // for the fracture model
  computeHistoryVariable(F_pos, F_neg);
  
  // residual of undamaged stress used for plasticity NR algorithm
  _residual_tensor = _pk2[_qp] - pk2_new;
}

// Jacobian for the Newton-Raphson crystal plasticity algorithm
// therefore, it does not include damage
void
CrystalPlasticityUndamagedStress::calculateJacobian()
{
  // may not need to cache the dfpinvdpk2 here. need to double check
  RankFourTensor dfedfpinv, deedfe, dfpinvdpk2, dfpinvdpk2_per_model;

  RankTwoTensor ffeiginv = _temporary_deformation_gradient * _inverse_eigenstrain_deformation_grad;

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
        dfedfpinv(i, j, k, j) = ffeiginv(i, k);

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _dislocation_models[i]->calculateTotalPlasticDeformationGradientDerivative(
        dfpinvdpk2_per_model,
        _inverse_plastic_deformation_grad_old,
        _inverse_eigenstrain_deformation_grad,
        _num_eigenstrains);
    dfpinvdpk2 += dfpinvdpk2_per_model;
  }
  
  _jacobian = RankFourTensor::IdentityFour() - (_elasticity_tensor[_qp] * deedfe * dfedfpinv * dfpinvdpk2);
}

// pk2_new is the undamaged stress
void
CrystalPlasticityUndamagedStress::computeStrainVolumetric(Real & F_pos, Real & F_neg, 
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
  
  // Undamaged elastic energy tensor and thermal energy
  RankTwoTensor undamaged_elastic_energy;
  RankTwoTensor undamaged_thermal_energy;
  
  // Positive and negative part of the Piola-Kirchhoff stress
  RankTwoTensor pk2_pos;
  RankTwoTensor pk2_neg;
  
  // inverse of the right Cauchy–Green deformation tensor (elastic)
  RankTwoTensor invce;
  
  // Thermal eigenstrain
  RankTwoTensor thermal_eigenstrain;
  
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
    a_pos_vol -= Kb * delta * _volumetric_thermal_expansion;
    a_neg_vol = 0.0;
	  
  } else { // compression
	  
    a_pos_vol = 0.0;
	a_neg_vol = 0.5 * Kb * delta * delta;
	a_neg_vol -= Kb * delta * _volumetric_thermal_expansion;
	
  }
  
  // Calculate isochoric free energy
  // Equation 17 in Grilli, Koslowski, 2019

  undamaged_elastic_energy = _elasticity_tensor[_qp] * ee;
  undamaged_elastic_energy = 0.5 * ee * undamaged_elastic_energy;
  
  undamaged_thermal_energy = _elasticity_tensor[_qp] * _thermal_eigenstrain;
  undamaged_thermal_energy = ee * undamaged_thermal_energy;

  a_pos_cpl = undamaged_elastic_energy.trace();
  a_pos_cpl -= 0.5 * Kb * delta * delta;
  
  a_pos_cpl -= undamaged_thermal_energy.trace();
  a_pos_cpl += Kb * delta * _volumetric_thermal_expansion;

  // Calculate positive and negative parts of the Piola-Kirchhoff stress
  // Equation 18 in Grilli, Koslowski, 2019
  if (Je >= 1.0) { // expansion
  
    pk2_pos = _elasticity_tensor[_qp] * (ee - _thermal_eigenstrain);
    pk2_neg = 0.0;
	
  } else { // compression
  
    invce = ce.inverse();
	
	pk2_neg = Je23 * Kb * delta * invce;
	
    pk2_pos = (-1.0) * pk2_neg;	
    pk2_pos += _elasticity_tensor[_qp] * (ee - _thermal_eigenstrain);

  }

  // Positive part of the stress is degraded by _D[_qp]
  _pk2_damaged[_qp] = (_D[_qp] * pk2_pos) + pk2_neg;
  pk2_new = pk2_pos + pk2_neg;
  
  // Positive and negative parts of the free energy
  // Equations 13 and 14 in Grilli, Koslowski, 2019
  // additionally, plastic work is included for damage
  F_pos = a_pos_vol + a_pos_cpl + _plastic_damage_prefactor * _plastic_work[_qp];
  F_neg = a_neg_vol;
  
  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  // This is used for stress equilibrium, therefore it has to account for the damage
  _dstress_dc[_qp] = pk2_pos * _dDdc[_qp];
  
  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history variable
  // This is used for phase field jacobian, therefore it has to account for the damage
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = pk2_pos * _dDdc[_qp];
  
}

// compute history variable and assign to _E
// which is used by the fracture model for damage growth
// Damage grows only because of the positive part of the elastic energy F_pos
void
CrystalPlasticityUndamagedStress::computeHistoryVariable(Real & F_pos, Real & F_neg)
{
  // Assign history variable
  Real hist_variable = _H_old[_qp];
  
  if (_use_snes_vi_solver) {
	  
    _H[_qp] = F_pos;
	  
  } else {
	  
    if (F_pos > _H_old[_qp])
      _H[_qp] = F_pos;
    else
      _H[_qp] = _H_old[_qp];
  }

  if (_use_current_hist)
    hist_variable = _H[_qp];

  // Elastic free energy density and derivatives
  // These are used by the phase field model, therefore they must include damage
  _E[_qp] = hist_variable * _D[_qp] + F_neg;
  _dEdc[_qp] = hist_variable * _dDdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp];
}

// update jacobian_mult by taking into account of the exact elasto-plastic tangent moduli
// it includes damage
void
CrystalPlasticityUndamagedStress::elastoPlasticTangentModuli(RankFourTensor & jacobian_mult)
{
  RankFourTensor tan_mod;
  RankTwoTensor pk2fet, fepk2, feiginvfpinv;
  RankFourTensor deedfe, dsigdpk2dfe, dfedf;
  
  const auto je = _elastic_deformation_gradient.det();

  // Fill in the matrix stiffness material property
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  usingTensorIndices(i_, j_, k_, l_);
  
  if (je >= 1.0) { // expansion
	  
    dsigdpk2dfe = _elastic_deformation_gradient.times<i_, k_, j_, l_>(_elastic_deformation_gradient) *
                  _D[_qp] * _elasticity_tensor[_qp] * deedfe;

  } else { // compression
	
    dsigdpk2dfe = _elastic_deformation_gradient.times<i_, k_, j_, l_>(_elastic_deformation_gradient) *
                  _elasticity_tensor[_qp] * deedfe;

  }
  
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

  if (je > 0.0)
    tan_mod /= je;

  feiginvfpinv = _inverse_eigenstrain_deformation_grad * _inverse_plastic_deformation_grad;
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
        dfedf(i, j, i, l) = feiginvfpinv(l, j);

  jacobian_mult = tan_mod * dfedf;
}
