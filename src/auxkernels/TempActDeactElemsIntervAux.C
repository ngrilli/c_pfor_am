// Nicolò Grilli
// 11 Gennaio 2021
// National University of Singapore

#include "TempActDeactElemsIntervAux.h"

#include <cmath>

registerMooseObject("TensorMechanicsApp", TempActDeactElemsIntervAux);

InputParameters
TempActDeactElemsIntervAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Read temperature due to laser scan during SLM from CFD simulation");
  params.addParam<std::string>("base_name", "Mechanical property base name");
  params.addParam<UserObjectName>("temperature_read_user_object",
                                  "The LaserTempReadFile "
                                  "GeneralUserObject to read element "
                                  "specific temperature values from file");
  params.addRequiredParam<Real>("temperature_time_step","Time interval between two temperature data field");
  params.addParam<Real>("melting_temperature_high", 1673.15, "Melting temperature (liquidus) = zero stiffness.");  
  params.addParam<Real>("melting_temperature_low", 1648.15, "Solidus = full stiffness.");
  params.addParam<Real>("gas_temperature_high", 298.1, "Lowest possible solid temperature = full stiffness.");
  params.addParam<Real>("gas_temperature_low", 298.0, "Gas temperature = zero stiffness.");
  params.addParam<Real>("deact_interval", 0.20, "Fraction of temperature_time_step for element deactivation");
  return params;
}

TempActDeactElemsIntervAux::TempActDeactElemsIntervAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
	_temperature_read_user_object(isParamValid("temperature_read_user_object")
                                  ? &getUserObject<LaserTempReadFile>("temperature_read_user_object")
                                  : nullptr),
    _temperature_time_step(getParam<Real>("temperature_time_step")),
	_melting_temperature_high(getParam<Real>("melting_temperature_high")),
	_melting_temperature_low(getParam<Real>("melting_temperature_low")),
	_gas_temperature_high(getParam<Real>("gas_temperature_high")),
	_gas_temperature_low(getParam<Real>("gas_temperature_low")),
	_deact_interval(getParam<Real>("deact_interval"))
{
}

// Calculate temperature read from file based on element index
// It can be used with GeneratedMesh
Real
TempActDeactElemsIntervAux::computeValue()
{
  Real TempValue = 0.0;
  Real TempValueNext = 0.0;
  unsigned int temperature_step;
  Real FracTimeStep = 0.0; // fraction of temperature time step completed, between 0 and 1
  
  // determine time step to be used from the CFD simulations
  temperature_step = floor(_t / _temperature_time_step);
  FracTimeStep = _t / _temperature_time_step - temperature_step;
	
  if (_temperature_read_user_object)
  {
    TempValue = _temperature_read_user_object->getData(_current_elem, temperature_step);
	TempValueNext = _temperature_read_user_object->getData(_current_elem, temperature_step+1);
	
	if (TempValue < _gas_temperature_high) { // from gas ...
	
		if (TempValueNext < _gas_temperature_high) { // ... to gas
		
			TempValue = _gas_temperature_low;
		
		} else if (TempValueNext <= _melting_temperature_low) { // ... to solid
		
			if (FracTimeStep > 1.0 - _deact_interval) {
				
				TempValue = ((1.0-FracTimeStep)/_deact_interval) * TempValue 
				          + ((_deact_interval-1.0+FracTimeStep)/_deact_interval) * TempValueNext;
				
			} // else TempValue remains the one read from file
			
		} else { // ... to liquid
		
			TempValue = _melting_temperature_high;
			
		}	

	} else if (TempValue <= _melting_temperature_low) { // from solid ...
		
		if (FracTimeStep > 1.0 - _deact_interval) {
		
			TempValue = ((1.0-FracTimeStep)/_deact_interval) * TempValue 
				      + ((_deact_interval-1.0+FracTimeStep)/_deact_interval) * TempValueNext;
			
		} // else TempValue remains the one read from file
		
	} else { // from liquid ...
	
		if (TempValueNext < _gas_temperature_high) { // ... to gas
		
			TempValue = _gas_temperature_low;
		
		} else if (TempValueNext <= _melting_temperature_low) { // ... to solid
		
			if (FracTimeStep > 1.0 - _deact_interval) {
			
				TempValue = ((1.0-FracTimeStep)/_deact_interval) * TempValue 
						  + ((_deact_interval-1.0+FracTimeStep)/_deact_interval) * TempValueNext;
			
			} // else TempValue remains the one read from file
			
		} else { // ... to liquid
			
			TempValue = _melting_temperature_high;
			
		}
	}
  }
  else
  {
	mooseError("Error in reading temperature file");
  }
  
  return TempValue;
}
