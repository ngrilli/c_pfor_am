// Nicolò Grilli
// Zhuohao Song
// Università di Bristol

#include "FunctionPathConicalHeatSource.h"

#include "Function.h"

registerMooseObject("c_pfor_amApp", FunctionPathConicalHeatSource);

InputParameters
FunctionPathConicalHeatSource::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<Real>("power_ellipsoid_1", "power of first ellipsoid");
  params.addRequiredParam<Real>("power_ellipsoid_2", "power of second ellipsoid");
  params.addRequiredParam<Real>("power_conical", "power of conical");
  params.addParam<Real>("efficiency_ellipsoid_1", 1, "process efficiency of first ellipsoid");
  params.addParam<Real>("efficiency_ellipsoid_2", 1, "process efficiency of second ellipsoid");
  params.addParam<Real>("efficiency_conical", 1, "process efficiency of conical");
  params.addRequiredParam<Real>("rx_ellipsoid_1", "effective transverse ellipsoid radius of first ellipsoid");
  params.addRequiredParam<Real>("ry_ellipsoid_1", "effective longitudinal ellipsoid radius of first ellipsoid");
  params.addRequiredParam<Real>("rz_ellipsoid_1", "effective depth ellipsoid radius of first ellipsoid");
  params.addRequiredParam<Real>("rx_ellipsoid_2", "effective transverse ellipsoid radius of second ellipsoid");
  params.addRequiredParam<Real>("ry_ellipsoid_2", "effective longitudinal ellipsoid radius of second ellipsoid");
  params.addRequiredParam<Real>("rz_ellipsoid_2", "effective depth ellipsoid radius of second ellipsoid");
  params.addRequiredParam<Real>("rx_conical", "effective transverse conical radius");
  params.addRequiredParam<Real>("ry_conical", "effective longitudinal conical radius");
  params.addRequiredParam<Real>("rz_conical", "effective depth conical radius: this should be the same as the layer thickness");
  params.addRequiredParam<Real>("cone_radius_decrease", "dimensionless quantity which corresponds to the axial decrease in the radius of the cone");
  params.addParam<Real>(
      "factor_ellipsoid_1", 1, "scaling factor that is multiplied to the heat source to adjust the intensity");
  params.addParam<Real>(
      "factor_ellipsoid_2", 1, "scaling factor that is multiplied to the heat source to adjust the intensity");
  params.addParam<Real>(
      "factor_conical", 1, "scaling factor that is multiplied to the heat source to adjust the intensity");
  params.addParam<FunctionName>(
      "function_x_ellipsoid_1", "0", "The x component of the center of the first ellipsoidal heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_y_ellipsoid_1", "0", "The y component of the center of the first ellipsoidal heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_z_ellipsoid_1", "0", "The z component of the center of the first ellipsoidal heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_x_ellipsoid_2", "0", "The x component of the center of the second ellipsoidal heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_y_ellipsoid_2", "0", "The y component of the center of the second ellipsoidal heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_z_ellipsoid_2", "0", "The z component of the center of the second ellipsoidal heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_x_conical", "0", "The x component of the center of the conical heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_y_conical", "0", "The y component of the center of the conical heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_z_conical", "0", "The z component of the center of the conical heating spot as a function of time: "
                                 "this should be the same as the back coordinate of the layer.");
  params.addClassDescription("Sum between two ellipsoidal heat source distribution and a conical heat source distribution.");

  return params;
}

FunctionPathConicalHeatSource::FunctionPathConicalHeatSource(const InputParameters & parameters)
  : Material(parameters),
    _P_ellipsoid_1(getParam<Real>("power_ellipsoid_1")),
    _P_ellipsoid_2(getParam<Real>("power_ellipsoid_2")),
    _P_conical(getParam<Real>("power_conical")),
    _eta_ellipsoid_1(getParam<Real>("efficiency_ellipsoid_1")),
    _eta_ellipsoid_2(getParam<Real>("efficiency_ellipsoid_2")),
    _eta_conical(getParam<Real>("efficiency_conical")),
    _rx_ellipsoid_1(getParam<Real>("rx_ellipsoid_1")),
    _ry_ellipsoid_1(getParam<Real>("ry_ellipsoid_1")),
    _rz_ellipsoid_1(getParam<Real>("rz_ellipsoid_1")),
    _rx_ellipsoid_2(getParam<Real>("rx_ellipsoid_2")),
    _ry_ellipsoid_2(getParam<Real>("ry_ellipsoid_2")),
    _rz_ellipsoid_2(getParam<Real>("rz_ellipsoid_2")),
    _rx_conical(getParam<Real>("rx_conical")),
    _ry_conical(getParam<Real>("ry_conical")),
    _rz_conical(getParam<Real>("rz_conical")),
    _cone_radius_decrease(getParam<Real>("cone_radius_decrease")),
    _f_ellipsoid_1(getParam<Real>("factor_ellipsoid_1")),
    _f_ellipsoid_2(getParam<Real>("factor_ellipsoid_2")),
    _f_conical(getParam<Real>("factor_conical")),
    _function_x_ellipsoid_1(getFunction("function_x_ellipsoid_1")),
    _function_y_ellipsoid_1(getFunction("function_y_ellipsoid_1")),
    _function_z_ellipsoid_1(getFunction("function_z_ellipsoid_1")),
    _function_x_ellipsoid_2(getFunction("function_x_ellipsoid_2")),
    _function_y_ellipsoid_2(getFunction("function_y_ellipsoid_2")),
    _function_z_ellipsoid_2(getFunction("function_z_ellipsoid_2")),
    _function_x_conical(getFunction("function_x_conical")),
    _function_y_conical(getFunction("function_y_conical")),
    _function_z_conical(getFunction("function_z_conical")),
    _volumetric_heat(declareADProperty<Real>("volumetric_heat"))
{
}

void
FunctionPathConicalHeatSource::computeQpProperties()
{
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);

  // center of the heat source: ellipsoid 1
  Real x_t_e1 = _function_x_ellipsoid_1.value(_t);
  Real y_t_e1 = _function_y_ellipsoid_1.value(_t);
  Real z_t_e1 = _function_z_ellipsoid_1.value(_t);

  // center of the heat source: ellipsoid 2
  Real x_t_e2 = _function_x_ellipsoid_2.value(_t);
  Real y_t_e2 = _function_y_ellipsoid_2.value(_t);
  Real z_t_e2 = _function_z_ellipsoid_2.value(_t);

  // center of the heat source: conical
  Real x_t_c = _function_x_conical.value(_t);
  Real y_t_c = _function_y_conical.value(_t);
  Real z_t_c = _function_z_conical.value(_t);

  // power of the heat sources
  Real q_e1;
  Real q_e2;
  Real q_c;

  q_e1 = 3.0 * std::sqrt(3.0) * _P_ellipsoid_1 * _eta_ellipsoid_1 * _f_ellipsoid_1 / (_rx_ellipsoid_1 * _ry_ellipsoid_1 * _rz_ellipsoid_1 * std::pow(libMesh::pi, 1.5)) *
         std::exp(-(3.0 * std::pow(x - x_t_e1, 2.0) / std::pow(_rx_ellipsoid_1, 2.0) +
                    3.0 * std::pow(y - y_t_e1, 2.0) / std::pow(_ry_ellipsoid_1, 2.0) +
                    3.0 * std::pow(z - z_t_e1, 2.0) / std::pow(_rz_ellipsoid_1, 2.0)));

  q_e2 = 3.0 * std::sqrt(3.0) * _P_ellipsoid_2 * _eta_ellipsoid_2 * _f_ellipsoid_2 / (_rx_ellipsoid_2 * _ry_ellipsoid_2 * _rz_ellipsoid_2 * std::pow(libMesh::pi, 1.5)) *
         std::exp(-(3.0 * std::pow(x - x_t_e2, 2.0) / std::pow(_rx_ellipsoid_2, 2.0) +
                    3.0 * std::pow(y - y_t_e2, 2.0) / std::pow(_ry_ellipsoid_2, 2.0) +
                    3.0 * std::pow(z - z_t_e2, 2.0) / std::pow(_rz_ellipsoid_2, 2.0)));

  q_c = 6.0 * _P_conical * _eta_conical * _f_conical / (_rx_conical * _ry_conical * _rz_conical * (1.0 + _cone_radius_decrease) * libMesh::pi) *
        std::exp(-(3.0 * std::pow(x - x_t_c, 2.0) / std::pow(_rx_conical, 2.0) +
                   3.0 * std::pow(y - y_t_c, 2.0) / std::pow(_ry_conical, 2.0)))
        * (1.0 - (1.0 - _cone_radius_decrease) * std::abs(z - z_t_c) / _rz_conical);

  _volumetric_heat[_qp] = q_e1 + q_e2 + q_c;
}
