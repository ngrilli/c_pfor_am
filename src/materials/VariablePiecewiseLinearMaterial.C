// Nicolò Grilli
// Università di Bristol
// 18 Aprile 2026

#include "VariablePiecewiseLinearMaterial.h"

registerMooseObject("c_pfor_amApp", VariablePiecewiseLinearMaterial);

InputParameters
VariablePiecewiseLinearMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "A material property that uses a piecewise linear function of a coupled variable");
  params.addRequiredCoupledVar("variable", "Variable used as input for the function");
  params.addRequiredParam<std::vector<Real>>("x", "The abscissa values");
  params.addRequiredParam<std::vector<Real>>("y", "The ordinate values");
  params.addRequiredParam<MaterialPropertyName>("property", "Output property");

  return params;
}

VariablePiecewiseLinearMaterial::VariablePiecewiseLinearMaterial(const InputParameters & parameters)
  : Material(parameters),
    _var(coupledValue("variable")),
    _x(getParam<std::vector<Real>>("x")),
    _y(getParam<std::vector<Real>>("y")),
    _f(declareProperty<Real>(getParam<MaterialPropertyName>("property")))
{
  if (_x.size() != _y.size())
    mooseError("x and y must have same size");
}

void
VariablePiecewiseLinearMaterial::computeQpProperties()
{
  const Real v = _var[_qp];  // value of the variable

  // clamp outside range
  if (v <= _x.front())
  {
    _f[_qp] = _y.front();
    return;
  }

  if (v >= _x.back())
  {
    _f[_qp] = _y.back();
    return;
  }

  // find interval
  for (unsigned int i = 0; i < _x.size() - 1; ++i)
  {
    if (v >= _x[i] && v < _x[i + 1])
    {
      const Real slope = (_y[i + 1] - _y[i]) / (_x[i + 1] - _x[i]);
      _f[_qp] = _y[i] + slope * (v - _x[i]);
      return;
    }
  }

  mooseWarning("This should never happen, but value is outside of x range and was not caught by clamping logic. Returning last y value as fallback.");
  _f[_qp] = _y.back(); // fallback (should never hit)
}
