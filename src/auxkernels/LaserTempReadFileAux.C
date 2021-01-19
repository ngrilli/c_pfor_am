// Nicolò Grilli
// 8 Novembre 2020
// National University of Singapore

#include "LaserTempReadFileAux.h"

#include <cmath>

registerMooseObject("TensorMechanicsApp", LaserTempReadFileAux);

InputParameters
LaserTempReadFileAux::validParams()
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
  params.addParam<Real>("reference_temperature",303.0,"reference temperature for thermal expansion");
  params.addParam<bool>("degrade_eigenstrain",false,"If liquid or gas, output room temperature to degrade eigenstrain");
  return params;
}

LaserTempReadFileAux::LaserTempReadFileAux(const InputParameters & parameters)
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
    _reference_temperature(getParam<Real>("reference_temperature")),
	_degrade_eigenstrain(getParam<bool>("degrade_eigenstrain"))
{
}

// Calculate temperature read from file based on element index
// It can be used with GeneratedMesh
Real
LaserTempReadFileAux::computeValue()
{
  Real TempValue = 0.0;
  Real TempValueNext = 0.0;
  Real TempValuePrevious = 0.0;
  unsigned int temperature_step;
  Real FracTimeStep = 0.0; // fraction of temperature time step completed, between 0 and 1
  
  unsigned int isSolidPrevious = 1; // to check if element has activated at current CFD step
  
  // determine time step to be used from the CFD simulations
  temperature_step = floor(_t / _temperature_time_step);
  FracTimeStep = _t / _temperature_time_step - temperature_step;
	
  if (_temperature_read_user_object)
  {
    TempValue = _temperature_read_user_object->getData(_current_elem, temperature_step);
	TempValueNext = _temperature_read_user_object->getData(_current_elem, temperature_step+1);

	// Limit temperature in the interval gas to liquidus
    // to avoid problem with the temperature dependencies
    // of elastic constants, CRSS, CTE
	
    TempValue = std::min(_melting_temperature_high,TempValue);
	TempValue = std::max(_gas_temperature_low,TempValue);
	
	TempValueNext = std::min(_melting_temperature_high,TempValueNext);
	TempValueNext = std::max(_gas_temperature_low,TempValueNext);
	
	if (_degrade_eigenstrain) {
		
	  if (temperature_step > 0) {
		  
		TempValuePrevious = _temperature_read_user_object->getData(_current_elem, temperature_step-1);
		
		TempValuePrevious = std::min(_melting_temperature_high,TempValuePrevious);
		TempValuePrevious = std::max(_gas_temperature_low,TempValuePrevious);
		
	    // check phase at the previous time step
	    if (TempValuePrevious < _gas_temperature_high) {
	      isSolidPrevious = 0;
	    } else if (TempValuePrevious <= _melting_temperature_low) {
	      isSolidPrevious = 1;
	    } else {
	      isSolidPrevious = 0;	
	    }	
		
	  }
	  
	  if (isSolidPrevious == 1) { // This element has not activated at this CFD step
		  
	  // If phase transition at the current or next CFD step
      // then transform linearly the temperature to room temperature
      // to degrade the eigenstrain to 0 in FiniteStrainCrystalPlasticityThermal

	  if (TempValue < _gas_temperature_high) { // from gas ...
	
		if (TempValueNext < _gas_temperature_high) { // ... to gas
		
			TempValue = _reference_temperature;
		
		} else if (TempValueNext <= _melting_temperature_low) { // ... to solid
		
		  TempValue = (1.0-FracTimeStep) * _reference_temperature 
				    + FracTimeStep * TempValueNext;
							
		} else { // ... to liquid
		
			TempValue = _reference_temperature;
			
		}	

	  } else if (TempValue <= _melting_temperature_low) { // from solid ...
	
	    if (TempValueNext < _gas_temperature_high) { // ... to gas
	
	      TempValue = (1.0-FracTimeStep) * TempValue 
		            + FracTimeStep * _reference_temperature;
					
		} else if (TempValueNext <= _melting_temperature_low) { // ... to solid
		  
 	      TempValue = (1.0-FracTimeStep) * TempValue 
		            + FracTimeStep * TempValueNext;         
		  
		} else { // ... to liquid
		
	      TempValue = (1.0-FracTimeStep) * TempValue 
		            + FracTimeStep * _reference_temperature;
					
		}
		
	  } else { // from liquid ...
	
		if (TempValueNext < _gas_temperature_high) { // ... to gas
		
		  TempValue = _reference_temperature;
		
		} else if (TempValueNext <= _melting_temperature_low) { // ... to solid
		
		  TempValue = (1.0-FracTimeStep) * _reference_temperature 
				    + FracTimeStep * TempValueNext;
			
		} else { // ... to liquid
			
		  TempValue = _reference_temperature;
			
		}
	  } 		  
		  
	  } else { // Element just activated: degrade thermal eigenstrain independently of current state
		
	    if (TempValueNext < _gas_temperature_high) { // ... to gas
		
		  TempValue = _reference_temperature;
		
		} else if (TempValueNext <= _melting_temperature_low) { // ... to solid
		
		  TempValue = (1.0-FracTimeStep) * _reference_temperature 
				    + FracTimeStep * TempValueNext;
			
		} else { // ... to liquid
		
		  TempValue = _reference_temperature;
			
		}
	  }
	} else {
	  // linear interpolation of the temperature in time
	  TempValue = (1.0 - FracTimeStep) * TempValue + FracTimeStep * TempValueNext;		
	}
  }
  else
  {
	mooseError("Error in reading temperature file");
  }
  
  return TempValue;
}
