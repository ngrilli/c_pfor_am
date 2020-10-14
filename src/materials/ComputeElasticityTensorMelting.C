// Nicolo Grilli
// ADSC Singapore
// 13 Ottobre 2020

#include "ComputeElasticityTensorMelting.h"
#include "RotationTensor.h"

registerMooseObject("TensorMechanicsApp", ComputeElasticityTensorMelting);

InputParameters
ComputeElasticityTensorMelting::validParams()
{
  InputParameters params = ComputeElasticityTensorCPGrain::validParams();
  params.addClassDescription("Compute an elasticity tensor for crystal plasticity. "
                             "Euler angles are read from Euler angles input file "
							 "and can be assigned to physical volumes in GMSH. "
							 "Melting is considered: stiffness is degraded when the "
							 "temperature increases from recrystallization to melting.");
  params.addParam<Real>("melting_temperature", 1300.0, "Melting temperature = zero stiffness.");  
  params.addParam<Real>("crystallisation_temperature", 800.0, "Crystallisation temperature = full stiffness.");
  return params;
}

ComputeElasticityTensorMelting::ComputeElasticityTensorMelting(const InputParameters & parameters)
  : ComputeElasticityTensorCPGrain(parameters),
	_melting_temperature(getParam<Real>("melting_temperature")),
	_crystallisation_temperature(getParam<Real>("crystallisation_temperature"))
{
	// _Cijkl is reinizialized to the unrotated state by the base class
}

void
ComputeElasticityTensorMelting::computeQpElasticityTensor()
{
  // Properties assigned at the beginning of every call to material calculation
  assignEulerAngles();

  _R.update(_Euler_angles_mat_prop[_qp]);

  _crysrot[_qp] = _R.transpose();
  
  // Apply temperature dependence on _Cijkl
  // and save results on _Temp_Cijkl
  temperatureDependence();
  
  // Apply degradation of the stiffness matrix
  // due to melting to _Temp_Cijkl
  // and save results on _Melt_Cijkl
  melting();
  
  _elasticity_tensor[_qp] = _Melt_Cijkl;

  _elasticity_tensor[_qp].rotate(_crysrot[_qp]);
}

void
ComputeElasticityTensorMelting::melting()
{	
  Real temp = _temp[_qp]; // Temperature
  Real melt_scale_factor;
  Real temp_interval; // Temperature interval between crystallisation and melting
  
  if (temp > _crystallisation_temperature) {
	  
	  // rescale the stiffness linearly with temperature
	  // in the interval [_crystallisation_temperature,_melting_temperature]
	  temp_interval = _melting_temperature - _crystallisation_temperature;
	  
	  if (temp_interval <= 0) {
		  mooseError("Melting temperature must be larger than crystallisation temperature");
	  }
	  
	  melt_scale_factor = std::min(1.0,(temp-_crystallisation_temperature)/temp_interval);
      _Melt_Cijkl = (1.0-melt_scale_factor) * _Temp_Cijkl;
	  
	  // Is a residual stiffness needed?
	  
  } else {
	  
	  _Melt_Cijkl = _Temp_Cijkl;
	  
  }
	
}


