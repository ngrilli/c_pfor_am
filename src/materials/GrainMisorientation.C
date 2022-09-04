// Nicolò Grilli
// University of Bristol
// 4 Settembre 2022

#include "GrainMisorientation.h"

registerMooseObject("TensorMechanicsApp", GrainMisorientation);

InputParameters
GrainMisorientation::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculate misorientation between neighbouring elements "
                             "of a structured mesh. ");
  return params;
}

GrainMisorientation::GrainMisorientation(const InputParameters & parameters)
  : Material(parameters),
//    _read_prop_user_object(isParamValid("read_prop_user_object")
//                               ? &getUserObject<PropertyReadFile>("read_prop_user_object")
//                               : nullptr),
    _misorientation(declareProperty<Real>("misorientation"))
{
}

void
GrainMisorientation::initQpStatefulProperties()
{
//  if (!_read_prop_user_object)
//  {
//    mooseError("read_prop_user_object not provided in GrainMisorientation");
//  }

  // TO DO: need to find a way to cycle over elements
  //_Euler_angles_mat_prop[_qp](0) = _read_prop_user_object->getData(_current_elem, 0);
  //_Euler_angles_mat_prop[_qp](1) = _read_prop_user_object->getData(_current_elem, 1);
  //_Euler_angles_mat_prop[_qp](2) = _read_prop_user_object->getData(_current_elem, 2);
  
  // This gives the number of elements
  // nobjects = _mesh.nElem();
  
  // This should provide the element pointer given the element id
  // auto el = _mesh.getMesh().query_elem_ptr(id);

  // TO DO: test the element id on a structured mesh 
  std::cout << this->_current_elem->id() << std::endl;
  
}

// misorientation is calculated only at the 
// first time step
void
GrainMisorientation::computeQpProperties()
{
}
