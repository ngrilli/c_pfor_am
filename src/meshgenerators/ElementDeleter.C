//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementDeleter.h"
#include "MooseMesh.h"

#include "libmesh/elem.h"

registerMooseObject("MooseApp", ElementDeleter);

template <>
InputParameters
validParams<ElementDeleter>()
{
    InputParameters params = validParams<ElementDeletionGeneratorBase>();
    params.addClassDescription(
        "Mesh modifier which removes elements with the element ID");
    params.addRequiredParam<std::vector<std::string>>("ids", "The IDs of element to be removed");
    return params;
}

ElementDeleter::ElementDeleter(const InputParameters& parameters)
  : ElementDeletionGeneratorBase(parameters), _ids((getParam<std::vector<std::string>>("ids")))
{
}

bool
ElementDeleter::shouldDelete(const Elem * elem)
{
  for (auto & charID : _ids){
    int delID = std::stoi(charID);
    if (elem->id() == delID){
      std::cout << elem->id() << std::endl;
      std::cout << "Element Deleted!" << std::endl;
      return true;
    }
  }
  auto elemID = elem->id();
  
  return false;
  //if (elemID == _id)
  //  std::cout << "ELEMENT DELETED" << std::endl;
  //return elem->id() == _id;
  
}

