// Nicolò Grilli
// Università di Bristol
// 29 Agosto 2024

#include "ReadFileIC.h"

registerMooseObject("c_pfor_amApp", ReadFileIC);

InputParameters
ReadFileIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addClassDescription("An initial condition that is read from file "
                             "and assigned to the nodes of a structured mesh. ");
  params.addParam<UserObjectName>("read_prop_user_object",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element "
                                  "specific property values from file. ");
  params.addParam<Real>("element_size", 1.0, "Element size in the structured mesh. ");
  return params;
}

ReadFileIC::ReadFileIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _read_prop_user_object(isParamValid("read_prop_user_object")
                               ? &getUserObject<PropertyReadFile>("read_prop_user_object")
                               : nullptr),
  _element_size(getParam<Real>("element_size"))
{
}

Real
ReadFileIC::value(const Point & p)
{
  // p(0) p(1) p(2) are the three coordinates of the point
  return p(0); //_func.value(_t, p);
}

RealGradient
ReadFileIC::gradient(const Point & p)
{
  return 0.0;
}
