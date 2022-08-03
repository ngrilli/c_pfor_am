// Nicolò Grilli
// University of Bristol
// 3 Agosto 2022

#include "TimeDepEllipsoidHeatSource.h"

#include "Function.h"

registerMooseObject("HeatConductionApp", TimeDepEllipsoidHeatSource);

InputParameters
TimeDepEllipsoidHeatSource::validParams()
{
  InputParameters params = FunctionPathEllipsoidHeatSource::validParams();
  params.addParam<FunctionName>(
      "function_t", "0", "The time function that is a prefactor of the space ellipsoid heat source.");
  params.addClassDescription("Double ellipsoid volumetric source heat with function path. "
                                             "A function of time is added for heat source ramp up and ramp down.");
  return params;
}

TimeDepEllipsoidHeatSource::TimeDepEllipsoidHeatSource(const InputParameters & parameters)
  : FunctionPathEllipsoidHeatSource(parameters),
    _function_t(getFunction("function_t"))
{
}

void
TimeDepEllipsoidHeatSource::computeQpProperties()
{
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);

  // center of the heat source
  Real x_t = _function_x.value(_t);
  Real y_t = _function_y.value(_t);
  Real z_t = _function_z.value(_t);
  
  // Time function for ramp up and ramp down
  Real ramp_up_t = _function_t.value(_t);

  _volumetric_heat[_qp] = 6.0 * ramp_up_t * std::sqrt(3.0) * _P * _eta * _f /
                          (_rx * _ry * _rz * std::pow(libMesh::pi, 1.5)) *
                          std::exp(-(3.0 * std::pow(x - x_t, 2.0) / std::pow(_rx, 2.0) +
                                     3.0 * std::pow(y - y_t, 2.0) / std::pow(_ry, 2.0) +
                                     3.0 * std::pow(z - z_t, 2.0) / std::pow(_rz, 2.0)));
}
