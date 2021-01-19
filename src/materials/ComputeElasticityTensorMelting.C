// Nicolo Grilli
// National University of Singapore
// 13 Novembre 2020

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
  params.addParam<Real>("residual_stiffness", 0.1, "Residual stiffness of gas and molten pool (percent).");
  params.addParam<UserObjectName>("temperature_read_user_object",
                                  "The LaserTempReadFile "
                                  "GeneralUserObject to read element "
                                  "specific temperature values from file");
  params.addParam<Real>("temperature_time_step",1.0,"Time interval between two temperature data field");
  params.addParam<bool>("activate_elems",false,"Using the element activation user object");
  return params;
}

ComputeElasticityTensorMelting::ComputeElasticityTensorMelting(const InputParameters & parameters)
  : ComputeElasticityTensorCPGrain(parameters),
	_melting_temperature_high(getParam<Real>("melting_temperature_high")),
	_melting_temperature_low(getParam<Real>("melting_temperature_low")),
	_gas_temperature_high(getParam<Real>("gas_temperature_high")),
	_gas_temperature_low(getParam<Real>("gas_temperature_low")),
	_residual_stiffness(getParam<Real>("residual_stiffness")),
	_temperature_read_user_object(isParamValid("temperature_read_user_object")
                                  ? &getUserObject<LaserTempReadFile>("temperature_read_user_object")
                                  : nullptr),
    _temperature_time_step(getParam<Real>("temperature_time_step")),
    _activate_elems(getParam<bool>("activate_elems"))
{
	// _Cijkl is reinizialized to the unrotated state by the base class
}

void
ComputeElasticityTensorMelting::computeQpElasticityTensor()
{
  Real temp = _temp[_qp]; // Temperature
  Real deltatemp;	
	
  // Properties assigned at the beginning of every call to material calculation
  assignEulerAngles();

  _R.update(_Euler_angles_mat_prop[_qp]);

  _crysrot[_qp] = _R.transpose();
  
  // Check phase at the current and next temperature time step
  checkPhase();
  
  if (_isSolidPrevious == 1 && _isSolid == 1 && _isSolidNext == 1) {
	  
	// Apply temperature dependence on _Cijkl
    // and save results on _Temp_Cijkl
    deltatemp = temp - _reference_temperature;
    temperatureDependence(deltatemp);
	
	_elasticity_tensor[_qp] = _Temp_Cijkl;
	
  } else {
	  
    // Apply degradation of the stiffness matrix
    // due to melting or gas to _Temp_Cijkl
    // and save results on _Melt_Cijkl
    melting(); 
	
    _elasticity_tensor[_qp] = _Melt_Cijkl;
  }
  
  _elasticity_tensor[_qp].rotate(_crysrot[_qp]);
}

void
ComputeElasticityTensorMelting::checkPhase()
{
  // determine time step to be used from the CFD simulations
  _temperature_step = floor(_t / _temperature_time_step);
  _FracTimeStep = _t / _temperature_time_step - _temperature_step;
  
  _isSolid = 0;
  _isLiquid = 0;
  _isGas = 0;
  _isSolidNext = 0;
  _isLiquidNext = 0;
  _isGasNext = 0;
  
  _isSolidPrevious = 1;
  
  if (_temperature_read_user_object)
  {
    _TempValue = _temperature_read_user_object->getData(_current_elem, _temperature_step);
	_TempValueNext = _temperature_read_user_object->getData(_current_elem, _temperature_step+1);
	
	if (_temperature_step > 0 && _activate_elems) { // check previous time step
      _TempValuePrevious = _temperature_read_user_object->getData(_current_elem, _temperature_step-1);
	  
	  _TempValuePrevious = std::min(_melting_temperature_high,_TempValuePrevious);
	  _TempValuePrevious = std::max(_gas_temperature_low,_TempValuePrevious);
	  
	  // check phase at the previous time step
	  if (_TempValuePrevious < _gas_temperature_high) {
	    _isSolidPrevious = 0;
	  } else if (_TempValuePrevious <= _melting_temperature_low) {
	    _isSolidPrevious = 1;
	  } else {
	    _isSolidPrevious = 0;	
	  }	
	}

	// Limit temperature in the interval gas to melting temperature
    // to avoid problem with the temperature dependencies
    // of elastic constants, CRSS, CTE
	
    _TempValue = std::min(_melting_temperature_high,_TempValue);
	_TempValue = std::max(_gas_temperature_low,_TempValue);
	
	_TempValueNext = std::min(_melting_temperature_high,_TempValueNext);
	_TempValueNext = std::max(_gas_temperature_low,_TempValueNext);
	
	// check phases, current and next, only one of the three flags will be activated
    if (_TempValue < _gas_temperature_high) {
	  _isGas = 1;
	} else if (_TempValue <= _melting_temperature_low) {
	  _isSolid = 1;
	} else {
	  _isLiquid = 1;	
	}	
	
	if (_TempValueNext < _gas_temperature_high) {
	  _isGasNext = 1;
	} else if (_TempValueNext <= _melting_temperature_low) {
	  _isSolidNext = 1;	  
	} else {
	  _isLiquidNext = 1;	
	}	
  } else {
	// Add code here to make this object working
    // when temperature is not read from external file
    // but it's just a variable	
	mooseError("Error in reading temperature file");	  
  }
}

void
ComputeElasticityTensorMelting::melting()
{	
  Real deltatemp;
  
  if (_isSolidPrevious == 1) { // case without element activation or previous solid CFD step
	  
    if (_isSolid == 1 && _isSolidNext == 0) { // becoming gas or liquid
  
  	  // start from temperature value at the last
	  // temperature time step
	  deltatemp = _TempValue - _reference_temperature;
	  temperatureDependence(deltatemp);
  
      _Melt_Cijkl = (1.0 - _FracTimeStep) * _Temp_Cijkl + _FracTimeStep * _residual_stiffness * _Cijkl;
		
    } else if (_isSolid == 0 && _isSolidNext == 1) { // back to solid
  
      // end at temperature value of the last
	  // temperature time step
	  deltatemp = _TempValueNext - _reference_temperature;
	  temperatureDependence(deltatemp);
  
      _Melt_Cijkl = (1.0 - _FracTimeStep) * _residual_stiffness * _Cijkl + _FracTimeStep * _Temp_Cijkl;
	  
    } else if (_isSolid == 0 && _isSolidNext == 0) { 
  
      // not solid at previous and next temperature time step:
	  // _Temp_Cijkl is never calculated in this case
	  _Melt_Cijkl = _residual_stiffness * _Cijkl;
	
    }	  
  } else { // case with element activation: stiffness degraded at current CFD step
  
      // In this case the phase of current CFD step is not important
	  if (_isSolidNext == 1) {
		  
        // end at temperature value of the last
	    // temperature time step
	    deltatemp = _TempValueNext - _reference_temperature;
	    temperatureDependence(deltatemp);
	    
        _Melt_Cijkl = (1.0 - _FracTimeStep) * _residual_stiffness * _Cijkl + _FracTimeStep * _Temp_Cijkl;
		
	  } else {
		  
        // not solid at previous and next temperature time step:
	    // _Temp_Cijkl is never calculated in this case
	    _Melt_Cijkl = _residual_stiffness * _Cijkl;
		  
	  }
  }
}


