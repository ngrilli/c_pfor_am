// Nicolò Grilli
// University of Bristol
// 4 Settembre 2022

#include "FileParsedMaterial.h"

registerMooseObject("TensorMechanicsApp", FileParsedMaterial);

InputParameters
FileParsedMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Parse a scalar material property from file "
                             "and assign to elements based on element number. ");
  params.addParam<std::string>("prop_name","arbitrary name of the material property. ");
  params.addParam<UserObjectName>("read_prop_user_object",
                                  "The ElementPropertyReadFile "
                                  "GeneralUserObject to read element "
                                  "specific property values from file");  
  return params;
}

FileParsedMaterial::FileParsedMaterial(const InputParameters & parameters)
  : Material(parameters),
    _prop_name(getParam<std::string>("prop_name")), // arbitrary name of the material property
    _property(declareProperty<Real>(_prop_name)), // Material property declaration
    _read_prop_user_object(isParamValid("read_prop_user_object")
                               ? &getUserObject<ElementPropertyReadFile>("read_prop_user_object")
                               : nullptr)
{
}

void
FileParsedMaterial::initQpStatefulProperties()
{
  computeQpProperties();
}

// The material property is only initialized and never updated
void
FileParsedMaterial::computeQpProperties()
{
  if (!_read_prop_user_object)
  {
    mooseError("read_prop_user_object not provided in FileParsedMaterial");
  }

  _property[_qp] = _read_prop_user_object->getData(_current_elem, 0);
}
