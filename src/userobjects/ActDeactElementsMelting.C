// Nicolo Grilli
// National University of Singapore
// 29 Dicembre 2020

#include "ActDeactElementsMelting.h"

registerMooseObject("MooseApp", ActDeactElementsMelting);

InputParameters
ActDeactElementsMelting::validParams()
{
  InputParameters params = ActDeactElementsUserObjectBase::validParams();
  
  params.addRequiredCoupledVar("temperature", "The temperature.");
  
  params.addParam<Real>("melting_temperature_low", 1648.15, "Delete element when liquid.");  
  params.addParam<Real>("gas_temperature_high", 298.1, "Delete element when gas.");

  return params;
}

ActDeactElementsMelting::ActDeactElementsMelting(const InputParameters & parameters)
  : ActDeactElementsUserObjectBase(parameters),
    _temperature(coupledValue("temperature")),
    _melting_temperature_low(
        declareRestartableData<Real>("melting_temperature_low", getParam<Real>("melting_temperature_low"))),
    _gas_temperature_high(
        declareRestartableData<Real>("gas_temperature_high", getParam<Real>("gas_temperature_high")))		
{
}

bool
ActDeactElementsMelting::isElementActivated()
{
  bool is_activated = true;
  Real avg_val = 0.0;

  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    avg_val += _temperature[qp];
  avg_val /= _qrule->n_points();
  
  if (avg_val > _melting_temperature_low) {
	  is_activated = false;
  }

  if (avg_val < _gas_temperature_high) {
	  is_activated = false;
  }

  return is_activated;
}
