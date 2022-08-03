// Nicolò Grilli
// University of Bristol
// 3 Agosto 2022

#pragma once

#include "FunctionPathEllipsoidHeatSource.h"

class Function;

/**
 * Double ellipsoid heat source distribution.
 * A function of time is added for heat source ramp up and ramp down.
 */
class TimeDepEllipsoidHeatSource : public FunctionPathEllipsoidHeatSource
{
public:
  static InputParameters validParams();

  TimeDepEllipsoidHeatSource(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// The time function that is a prefactor of the space ellipsoid heat source.
  const Function & _function_t;

};
