// Nicolò Grilli
// Parsa Esmati
// Zhuohao Song
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
  params.addParam<Real>("element_size_x", 1.0, "Element size along x in the structured mesh. ");
  params.addParam<Real>("element_size_y", 1.0, "Element size along y in the structured mesh. ");
  params.addParam<unsigned int>("nx", 1, "Number of nodes along the x axis in the structured mesh. ");
  params.addParam<unsigned int>("ny", 1, "Number of nodes along the y axis in the structured mesh. ");
  params.addParam<unsigned int>("op", 0, "Phase field index used to select the column in the file. ");
  params.addRequiredParam<unsigned int>("op_num", "Specifies the total number of phase fields to create. ");
  return params;
}

ReadFileIC::ReadFileIC(const InputParameters & parameters)
  : InitialCondition(parameters),
  _ic_file_name(getParam<FileName>("ic_file_name")),
  _element_size(getParam<Real>("element_size")),
  _different_xy_element_size(isParamValid("element_size_x") && isParamValid("element_size_y")),
  _element_size_x(getParam<Real>("element_size_x")),
  _element_size_y(getParam<Real>("element_size_y")),
  _nx(getParam<unsigned int>("nx")),
  _ny(getParam<unsigned int>("ny")),
  _op(getParam<unsigned int>("op")),
  _op_num(getParam<unsigned int>("op_num"))
{
  getFileData();
}

Real
ReadFileIC::value(const Point & p)
{
  unsigned int node_index;
  Real element_size_x;
  Real element_size_y;
  
  // check if mesh has different element sizes in x and y directions
  if (_different_xy_element_size) {
    element_size_x = _element_size_x;
    element_size_y = _element_size_y;
  } else {
    element_size_x = _element_size;
    element_size_y = _element_size;
  }
  
  // p(0) p(1) p(2) are the three coordinates of the point
  // Retrieve element index (0-based) from node coordinate  
  node_index = std::floor(p(0) / element_size_x) + _nx * std::floor(p(1) / element_size_y);
  
  // Retrieve the value of the variable from data structure
  return _IC_data[node_index][_op];
}

RealGradient
ReadFileIC::gradient(const Point & p)
{
  return 0.0;
}

void
ReadFileIC::getFileData()
{
  // read in the initial variable values from file
  MooseUtils::DelimitedFileReader _reader(_ic_file_name);
  _reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::ROWS);
  _reader.read();
  
  // check the size of the input corresponds to total number of elements
  if (_reader.getData().size() != _nx * _ny)
    paramError(
        "nx",
        "The number of rows in the ic_fileshould match the total number of elements. ");
        
  // resize data structure to store IC read from file
  _IC_data.resize(_nx * _ny);
  // and set the first index of each inner vector
  for (unsigned int i = 0; i < (_nx * _ny); ++i)
    _IC_data[i].resize(_op_num);
    
  for (const auto i : make_range(_nx * _ny))
  {
    for (const auto j : make_range(_op_num))
    {
      _IC_data[i][j] = _reader.getData(i)[j];
    }
  }
}
