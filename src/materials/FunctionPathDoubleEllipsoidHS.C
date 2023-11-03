// Nicolò Grilli
// Valentin Godard
// University of Bristol
// 2 Novembre 2023

#include "FunctionPathDoubleEllipsoidHS.h"

#include "Function.h"

registerMooseObject("HeatConductionApp", FunctionPathDoubleEllipsoidHS);

InputParameters
FunctionPathDoubleEllipsoidHS::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<Real>("power", "power");
  params.addParam<Real>("efficiency", 1, "process efficiency");
  params.addRequiredParam<Real>("a_f", "ellipsoid axis along the scanning direction: front");
  params.addRequiredParam<Real>("a_r", "ellipsoid axis along the scanning direction: rear");
  params.addRequiredParam<Real>("b", "ellipsoid axis perpendicular to the scanning direction");
  params.addRequiredParam<Real>("c", "ellipsoid axis perpendicular to the scanning direction");
  params.addParam<Real>(
      "factor_front", 1, "scaling factor that is multiplied to the heat source to adjust the intensity: front");
  params.addParam<Real>(
      "factor_rear", 1, "scaling factor that is multiplied to the heat source to adjust the intensity: rear");
  params.addParam<FunctionName>(
      "function_x", "0", "The x component of the center of the heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_y", "0", "The y component of the center of the heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_z", "0", "The z component of the center of the heating spot as a function of time");
  params.addClassDescription("Double ellipsoid volumetric source heat with function path based on: "
                             "O. Mokrov, M. Simon, A. Schiebahn, U. Reisgen, "
                             "A fine modification of the double ellipsoid heat source, "
                             "Mathematical Modelling of Weld Phenomena 12");

  return params;
}

FunctionPathDoubleEllipsoidHS::FunctionPathDoubleEllipsoidHS(const InputParameters & parameters)
  : Material(parameters),
    _P(getParam<Real>("power")),
    _eta(getParam<Real>("efficiency")),
    _a_f(getParam<Real>("a_f")),
    _a_r(getParam<Real>("a_r")),
    _b(getParam<Real>("b")),
    _c(getParam<Real>("c")),
    
    _f_f(getParam<Real>("factor_front")),
    _f_r(getParam<Real>("factor_rear")),
    
    _function_x(getFunction("function_x")),
    _function_y(getFunction("function_y")),
    _function_z(getFunction("function_z")),
    _volumetric_heat(declareADProperty<Real>("volumetric_heat"))
{
}

void
FunctionPathDoubleEllipsoidHS::computeQpProperties()
{
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);

  // center of the heat source
  Real x_t = _function_x.value(_t);
  Real y_t = _function_y.value(_t);
  Real z_t = _function_z.value(_t);
  
  // power of the front and rear heat sources
  Real q_f;
  Real q_r;

  q_f = 3.0 * std::sqrt(3.0) * _P * _eta * _f_f / (_b * _c * _a_f * std::pow(libMesh::pi, 1.5)) *
        std::exp(-3.0 * std::pow(x - x_t, 2.0) / std::pow(_a_f, 2.0) 
                 -3.0 * std::pow(y - y_t, 2.0) / std::pow(_b, 2.0) 
                 -3.0 * std::pow(z - z_t, 2.0) / std::pow(_c, 2.0));
  
  q_r = 3.0 * std::sqrt(3.0) * _P * _eta * _f_r / (_b * _c * _a_r * std::pow(libMesh::pi, 1.5)) *
        std::exp(-3.0 * std::pow(x - x_t, 2.0) / std::pow(_a_r, 2.0) 
                 -3.0 * std::pow(y - y_t, 2.0) / std::pow(_b,2.0) 
                 -3.0 * std::pow(z - z_t, 2.0) / std::pow(_c, 2.0));
                                     
  _volumetric_heat[_qp] = q_f + q_r;
}
