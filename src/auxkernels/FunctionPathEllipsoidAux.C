// Nicolò Grilli
// University of Bristol
// 7 Agosto 2022

#include "FunctionPathEllipsoidAux.h"

#include "Function.h"

registerMooseObject("MooseApp", FunctionPathEllipsoidAux);

InputParameters
FunctionPathEllipsoidAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("This AuxKernel increases an AuxVariable from 0 to 1 if the qp is inside "
                                             "an ellipsoid that is moving according to a path provided through a function. "
											 "It can be applied to a level set variable "
                                             "to simulate the material deposition during wire arc additive manufacturing "
											 "together with ActDeactElementsCoupled.");
  params.addRequiredCoupledVar("level_set_var", "The AuxVariable representing the level set.");
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

FunctionPathEllipsoidAux::FunctionPathEllipsoidAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _level_set_var(coupledValue("level_set_var")),
    _low_level_set_var(getParam<Real>("low_level_set_var")),
	_high_level_set_var(getParam<Real>("high_level_set_var")),
    _rx(getParam<Real>("rx")),
    _ry(getParam<Real>("ry")),
    _rz(getParam<Real>("rz")),
    _function_x(getFunction("function_x")),
    _function_y(getFunction("function_y")),
    _function_z(getFunction("function_z")),
	_level_set_activation_threshold(getParam<Real>("level_set_activation_threshold"))
{
}

Real
FunctionPathEllipsoidAux::computeValue()
{
  // value of the level set variable at the previous time step
  Real old_level_set = _u[_qp];
  
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);
  
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
	  
	  return _high_level_set_var;
	  
  } else {
	  
    if (old_level_set > _low_level_set_var) { // this was already activated
		 
      return _high_level_set_var;
		 
    } else { // this remains inactive

      return _low_level_set_var;

    }	  
	
  }
}
