// Nicolò Grilli
// Edward Horton
// University of Bristol
// 12 Aprile 2023

#include "ComputeCPStressR5CreepDamage.h"

registerMooseObject("TensorMechanicsApp", ComputeCPStressR5CreepDamage);

InputParameters
ComputeCPStressR5CreepDamage::validParams()
{
  InputParameters params = ComputeCrystalPlasticityStressDamage::validParams();
  params.addClassDescription(
      "This is a variant of ComputeCrystalPlasticityStressDamage in which the fracture energy "
      "is degraded based on the R5 creep damage criterion, which is based on "
      "ductility exhaustion theory."
	  "The damage model is reported in: "
	  "M. W. Spindler, "
	  "The prediction of creep damage in type 347 weld metal. "
	  "Part I: the determination of material properties from creep and tensile tests "
	  "International Journal of Pressure Vessels and Piping "
	  "Volume 82, Issue 3, March 2005, Pages 175-184 ");
  params.addParam<Real>("residual_creep_degradation", 1e-3, "Minimum residual creep degradation.");
  return params;
}

ComputeCPStressR5CreepDamage::ComputeCPStressR5CreepDamage(
    const InputParameters & parameters)
  : ComputeCrystalPlasticityStressDamage(parameters),
  _residual_creep_degradation(getParam<Real>("residual_creep_degradation")),
  _creep_degradation(declareProperty<Real>("creep_degradation")),
  _creep_degradation_old(getMaterialPropertyOld<Real>("creep_degradation"))
{
  _convergence_failed = false;
}

// Initialize internal state variables
void
ComputeCPStressR5CreepDamage::initQpStatefulProperties()
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
  
  // Initialize creep degradation function
  _creep_degradation[_qp] = 0.0;
}

// This is how a function from the dislocation model is called
// if some variable in the dislocation model are necessary
// to calculate the creep degradation
// _dislocation_models[i]->calculateFlowDirection(_crysrot[_qp]);

// Functions that update state variables 
// should be called here

// Fpos and Fneg should be calculated here for more efficiency

void
ComputeCPStressR5CreepDamage::postSolveQp(RankTwoTensor & cauchy_stress,
                                          RankFourTensor & jacobian_mult)
{
  cauchy_stress = _elastic_deformation_gradient * _pk2[_qp] *
                  _elastic_deformation_gradient.transpose() / _elastic_deformation_gradient.det();

  calcTangentModuli(jacobian_mult);
  
  // Calculate increment of Fp over _dt
  _fp_increment[_qp] = _plastic_deformation_gradient[_qp] - _plastic_deformation_gradient_old[_qp];
  
  // update plastic work
  // updatePlasticWork();

  _total_lagrangian_strain[_qp] =
      _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] -
      RankTwoTensor::Identity();
  _total_lagrangian_strain[_qp] = _total_lagrangian_strain[_qp] * 0.5;

  // Calculate crystal rotation to track separately
  RankTwoTensor rot;
  _elastic_deformation_gradient.getRUDecompositionRotation(rot);
  _updated_rotation[_qp] = rot * _crysrot[_qp];
}

// compute history variable and assign to _E
// which is used by the fracture model for damage growth
// Damage grows only because of the positive part of the elastic energy F_pos
void
ComputeCPStressR5CreepDamage::computeHistoryVariable(Real & F_pos, Real & F_neg)
{
  // Assign history variable
  Real hist_variable = _H_old[_qp];
  
  // _use_snes_vi_solver option not implemented

  if (F_pos > _H_old[_qp])
    _H[_qp] = F_pos;
  else
    _H[_qp] = _H_old[_qp];

  if (_use_current_hist)
    hist_variable = _H[_qp];

  // _barrier not implemented

  // Elastic free energy density and derivatives
  _E[_qp] = hist_variable * _D[_qp] + F_neg;
  _dEdc[_qp] = hist_variable * _dDdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp];

}
