// Nicolo Grilli
// National University of Singapore
// 10 Novembre 2020

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
							 "temperature increases above melting or below gas temperature.");
  params.addParam<Real>("melting_temperature_high", 1673.15, "Melting temperature (liquidus) = zero stiffness.");  
  params.addParam<Real>("melting_temperature_low", 1648.15, "Solidus = full stiffness.");
  params.addParam<Real>("gas_temperature_high", 298.1, "Lowest possible solid temperature = full stiffness.");
  params.addParam<Real>("gas_temperature_low", 298.0, "Gas temperature = zero stiffness.");
  return params;
}

ComputeElasticityTensorMelting::ComputeElasticityTensorMelting(const InputParameters & parameters)
  : ComputeElasticityTensorCPGrain(parameters),
	_melting_temperature_high(getParam<Real>("melting_temperature_high")),
	_melting_temperature_low(getParam<Real>("melting_temperature_low")),
	_gas_temperature_high(getParam<Real>("gas_temperature_high")),
	_gas_temperature_low(getParam<Real>("gas_temperature_low"))
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
  Real temp_interval; // Temperature interval between full and zero stiffness
  
  if (_gas_temperature_high > _melting_temperature_low) {
	  mooseError("Environment gas temperature cannot be higher than melting temperature");
  }
  
  if (temp > _melting_temperature_low) {
	  
	  // rescale the stiffness linearly with temperature
	  // in the interval [_melting_temperature_low,_melting_temperature_high]
	  temp_interval = _melting_temperature_high - _melting_temperature_low;
	  
	  if (temp_interval <= 0) {
		  mooseError("Liquidus temperature must be larger than solidus temperature");
	  }
	  
	  melt_scale_factor = std::min(1.0,(temp-_melting_temperature_low)/temp_interval);
      _Melt_Cijkl = (1.0-melt_scale_factor) * _Temp_Cijkl;
	  
	  // Is a residual stiffness needed?
	  
  } else if (temp < _gas_temperature_high) {
	  
	  // rescale the stiffness linearly with temperature
	  // in the interval [_gas_temperature_low,_gas_temperature_high]
	  temp_interval = _gas_temperature_high - _gas_temperature_low;
	  
	  if (temp_interval <= 0) {
		  mooseError("Gas temperatures low and high inverted");
	  }	  
	  
	  melt_scale_factor = std::max(0.0,(temp-_gas_temperature_low)/temp_interval);
	  _Melt_Cijkl = (1.0-melt_scale_factor) * _Temp_Cijkl;
	  
	  // Is a residual stiffness needed?
	  
  } else {
	  _Melt_Cijkl = _Temp_Cijkl;
  }
}


