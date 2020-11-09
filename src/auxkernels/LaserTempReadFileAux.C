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
  return params;
}

LaserTempReadFileAux::LaserTempReadFileAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
	_temperature_read_user_object(isParamValid("temperature_read_user_object")
                                  ? &getUserObject<LaserTempReadFile>("temperature_read_user_object")
                                  : nullptr),
    _temperature_time_step(getParam<Real>("temperature_time_step"))
{
}

// Calculate temperature read from file based on element index
// It can be used with GeneratedMesh
Real
LaserTempReadFileAux::computeValue()
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
	// linear interpolation of the temperature in time
	TempValue = (1.0 - FracTimeStep) * TempValue + FracTimeStep * TempValueNext;
  }
  else
  {
	mooseError("Error in reading temperature file");
  }
  
  return TempValue;
}
