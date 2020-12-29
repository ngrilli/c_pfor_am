// Nicolo Grilli
// National University of Singapore
// 29 Dicembre 2020

#pragma once

#include "ActDeactElementsUserObjectBase.h"

class ActDeactElementsMelting : public ActDeactElementsUserObjectBase
{
public:
  static InputParameters validParams();

  ActDeactElementsMelting(const InputParameters & parameters);

  virtual bool isElementActivated() override;

protected:
  /// temperature value to decide wether an element would be activated
  const VariableValue & _temperature;
  /// melting temperature above which elements are deactivated
  const Real _melting_temperature_low;
  /// gas temperature below which elements are deactivated
  const Real _gas_temperature_high;
};
