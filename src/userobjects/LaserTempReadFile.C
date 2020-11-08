// Nicolo Grilli
// National University of Singapore
// 8 Novembre 2020

#include "LaserTempReadFile.h"
#include "MooseRandom.h"
#include "MooseMesh.h"

#include <fstream>

registerMooseObject("TensorMechanicsApp", LaserTempReadFile);

InputParameters
LaserTempReadFile::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("User Object to read temperature data from an external file and assign "
                             "to elements.");
  params.addParam<FileName>("temperature_file_name","", "Name of the temperature file");
  params.addRequiredParam<unsigned int>("temperature_num_step","Number of temperature data field in time");
  return params;
}

LaserTempReadFile::LaserTempReadFile(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _temperature_file_name(getParam<FileName>("temperature_file_name")),
	_temperature_num_step(getParam<unsigned int>("temperature_num_step")),
    _mesh(_fe_problem.mesh())
{
  _nelem = _mesh.nElem();

  for (unsigned int i = 0; i < LIBMESH_DIM; i++)
  {
    _bottom_left(i) = _mesh.getMinInDimension(i);
    _top_right(i) = _mesh.getMaxInDimension(i);
    _range(i) = _top_right(i) - _bottom_left(i);
  }

  _max_range = _range(0);
  for (unsigned int i = 1; i < LIBMESH_DIM; i++)
    if (_range(i) > _max_range)
      _max_range = _range(i);

  readElementData();
}

void
LaserTempReadFile::readElementData()
{
  _data.resize(_nelem * _temperature_num_step);

  MooseUtils::checkFileReadable(_temperature_file_name);

  std::ifstream file_prop;
  file_prop.open(_temperature_file_name.c_str());

  for (unsigned int j = 0; j < _temperature_num_step; j++)
    for (unsigned int i = 0; i < _nelem; i++)
      if (!(file_prop >> _data[i + j * _nelem]))
        mooseError("Error LaserTempReadFile: Premature end of temperature file");

  file_prop.close();
}

Real
LaserTempReadFile::getData(const Elem * elem, unsigned int temperature_step) const
{
  unsigned int jelem = elem->id();
    
  mooseAssert(jelem < _nelem,
              "Error LaserTempReadFile: Element "
                  << jelem << " greater than total number of element in mesh " << _nelem);
  mooseAssert(temperature_step < _temperature_num_step,
              "Error LaserTempReadFile: Time step number "
                  << temperature_step << " greater than total number of temperature time steps " 
				  << _temperature_num_step);
  return _data[jelem + temperature_step * _nelem];
}
