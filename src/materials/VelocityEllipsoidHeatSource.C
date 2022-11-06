// Nicolò Grilli
// University of Bristol
// 6 Novembre 2022

#include "VelocityEllipsoidHeatSource.h"

#include "Function.h" // is this necessary?

registerMooseObject("HeatConductionApp", VelocityEllipsoidHeatSource);

InputParameters
VelocityEllipsoidHeatSource::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Double ellipsoid volumetric source heat, the motion is determined "
                             "by input velocity, starting positions and a postprocessor. ");
  params.addRequiredParam<Real>("power", "power");
  params.addParam<Real>("efficiency", 1, "process efficiency");
  params.addRequiredParam<Real>("rx", "effective transverse ellipsoid radius");
  params.addRequiredParam<Real>("ry", "effective longitudinal ellipsoid radius");
  params.addRequiredParam<Real>("rz", "effective depth ellipsoid radius");
  params.addParam<Real>(
      "factor", 1, "scaling factor that is multiplied to the heat source to adjust the intensity");
  params.addRequiredParam<RealVectorValue>("velocity", "Velocity vector");
  
  // Every time the postprocessor condition is satisfied, the heat source is moved to the next set of coordinates
  params.addRequiredParam<std::vector<Real>>("init_x_coords", "Initial values of x coordinates of the heat source");
  params.addRequiredParam<std::vector<Real>>("init_y_coords", "Initial values of y coordinates of the heat source");
  params.addRequiredParam<std::vector<Real>>("init_z_coords", "Initial values of z coordinates of the heat source");
  
  params.addRequiredParam<PostprocessorName>("temperature_pp","Postprocessor with temperature value to determine heat source motion.");
      
      
      
//  params.addParam<FunctionName>(
//      "function_x", "0", "The x component of the center of the heating spot as a function of time");
//  params.addParam<FunctionName>(
//      "function_y", "0", "The y component of the center of the heating spot as a function of time");
//  params.addParam<FunctionName>(
//      "function_z", "0", "The z component of the center of the heating spot as a function of time");


  return params;
}

VelocityEllipsoidHeatSource::VelocityEllipsoidHeatSource(const InputParameters & parameters)
  : Material(parameters),
    _P(getParam<Real>("power")),
    _eta(getParam<Real>("efficiency")),
    _rx(getParam<Real>("rx")),
    _ry(getParam<Real>("ry")),
    _rz(getParam<Real>("rz")),
    _f(getParam<Real>("factor")),
    _velocity(getParam<RealVectorValue>("velocity")), // Scanning speed vector
    
    // Initial values of the coordinates of the heat source
    _init_x_coords(getParam<std::vector<Real>>("init_x_coords")),
    _init_y_coords(getParam<std::vector<Real>>("init_y_coords")),
    _init_z_coords(getParam<std::vector<Real>>("init_z_coords")),
    
    // Postprocess with temperature value
    _temperature_pp(getPostprocessorValue("temperature_pp")),
    
//    _function_x(getFunction("function_x")),
//    _function_y(getFunction("function_y")),
//    _function_z(getFunction("function_z")),
    
    _volumetric_heat(declareADProperty<Real>("volumetric_heat"))
{
}

void
VelocityEllipsoidHeatSource::computeQpProperties()
{
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);

  // center of the heat source
  Real x_t = 0.0; //_function_x.value(_t);
  Real y_t = 0.0; // _function_y.value(_t);
  Real z_t = 0.0; // _function_z.value(_t);

  _volumetric_heat[_qp] = 6.0 * std::sqrt(3.0) * _P * _eta * _f /
                          (_rx * _ry * _rz * std::pow(libMesh::pi, 1.5)) *
                          std::exp(-(3.0 * std::pow(x - x_t, 2.0) / std::pow(_rx, 2.0) +
                                     3.0 * std::pow(y - y_t, 2.0) / std::pow(_ry, 2.0) +
                                     3.0 * std::pow(z - z_t, 2.0) / std::pow(_rz, 2.0)));
}
