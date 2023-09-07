// Nicolò Grilli
// Università di Bristol
// 3 Settembre 2023

#include "FunctionPathEllipsoidMaterial.h"

registerMooseObject("TensorMechanicsApp", FunctionPathEllipsoidMaterial);

InputParameters
FunctionPathEllipsoidMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "This class increases a material property from 0 to 1 if the element is inside an ellipsoid");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addParam<Real>("low_level_set_var", 0.0, "The lowest value of the level set variable.");
  params.addParam<Real>("high_level_set_var", 1.0, "The highest value of the level set variable.");
  params.addRequiredParam<Real>("rx", "effective transverse ellipsoid radius");
  params.addRequiredParam<Real>("ry", "effective longitudinal ellipsoid radius");
  params.addRequiredParam<Real>("rz", "effective depth ellipsoid radius");
  params.addParam<FunctionName>(
      "function_x", "0", "The x component of the center of the heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_y", "0", "The y component of the center of the heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_z", "0", "The z component of the center of the heating spot as a function of time");
  params.addParam<Real>("level_set_activation_threshold", 0.5, "Threshold value of the ellipsoid function "
                                                               "that activates the level set.");
  return params;
}

FunctionPathEllipsoidMaterial::FunctionPathEllipsoidMaterial(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _low_level_set_var(getParam<Real>("low_level_set_var")),
	_high_level_set_var(getParam<Real>("high_level_set_var")),
    _rx(getParam<Real>("rx")),
    _ry(getParam<Real>("ry")),
    _rz(getParam<Real>("rz")),
    _function_x(getFunction("function_x")),
    _function_y(getFunction("function_y")),
    _function_z(getFunction("function_z")),
	_level_set_activation_threshold(getParam<Real>("level_set_activation_threshold")),
    _level_set(declareProperty<Real>(_base_name + "level_set")),
    _level_set_old(getMaterialPropertyOld<Real>(_base_name + "level_set"))
{
}

void
FunctionPathEllipsoidMaterial::initQpStatefulProperties()
{
  _level_set[_qp] = 0.0;
}

void
FunctionPathEllipsoidMaterial::computeQpProperties()
{
  // value of the level set variable at the previous time step
  Real old_level_set = _level_set_old[_qp];
  
  const Real x = _q_point[_qp](0);
  const Real y = _q_point[_qp](1);
  const Real z = _q_point[_qp](2);
  
  // center of the ellipsoidal heat source
  Real x_t = _function_x.value(_t);
  Real y_t = _function_y.value(_t);
  Real z_t = _function_z.value(_t);
  
  // ellipsoid function value
  Real val;
  
  val = 6.0 * std::sqrt(3.0) /
          (_rx * _ry * _rz * std::pow(libMesh::pi, 1.5)) *
          std::exp(-(3.0 * std::pow(x - x_t, 2.0) / std::pow(_rx, 2.0) +
                          3.0 * std::pow(y - y_t, 2.0) / std::pow(_ry, 2.0) +
                           3.0 * std::pow(z - z_t, 2.0) / std::pow(_rz, 2.0)));

  if (val > _level_set_activation_threshold) { // heat source activating this _qp
	  
    _level_set[_qp] = _high_level_set_var;
	  
  } else {
	  
    if (old_level_set > _low_level_set_var) { // this was already activated
		 
      _level_set[_qp] = _high_level_set_var;
		 
    } else { // this remains inactive

      _level_set[_qp] = _low_level_set_var;

    }	  
	
  }
}
