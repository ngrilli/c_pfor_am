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
  params.addRequiredParam<FileName>(
      "ic_file_name",
      "Name of the file containing the initial conditions of a variable on a structured mesh. ");
  params.addParam<Real>("element_size", 1.0, "Element size in the structured mesh. ");
  params.addParam<unsigned int>("nx", 1, "Number of elements along the x axis in the structured mesh. ");
  params.addParam<unsigned int>("ny", 1, "Number of elements along the y axis in the structured mesh. ");
  return params;
}

ReadFileIC::ReadFileIC(const InputParameters & parameters)
  : InitialCondition(parameters),
  _ic_file_name(getParam<FileName>("ic_file_name")),
  _element_size(getParam<Real>("element_size")),
  _nx(getParam<unsigned int>("nx")),
  _ny(getParam<unsigned int>("ny"))
{
}

Real
ReadFileIC::value(const Point & p)
{
  unsigned int element_index;
  
  // p(0) p(1) p(2) are the three coordinates of the point
  // Retrieve element index (0-based) from node coordinate  
  element_index = std::floor(p(0) / _element_size) + _nx * std::floor(p(1) / _element_size);
  
  // read in the initial variable values from file
  MooseUtils::DelimitedFileReader _reader(_ic_file_name);
  _reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::ROWS);
  _reader.read();
  
  // check the size of the input corresponds to total number of elements
  if (_reader.getData().size() != _nx * _ny)
    paramError(
        "nx",
        "The number of rows in the ic_fileshould match the total number of elements. ");
  
  // Retrieve the value of the variable from file reader
  return _reader.getData(element_index)[0];
}

RealGradient
ReadFileIC::gradient(const Point & p)
{
  return 0.0;
}
