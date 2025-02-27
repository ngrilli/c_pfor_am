// Nicolò Grilli
// Università di Bristol
// Fernando Valiente Dies
// ANSTO
// 26 Febbraio 2025

#include "Eta2EulerAnglesMaterial.h"

registerMooseObject("c_pfor_amApp", Eta2EulerAnglesMaterial);

InputParameters
Eta2EulerAnglesMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Output euler angles from eta_i phase fields");
  params.addParam<FileName>(
      "Euler_angles_file_name","",
      "Name of the file containing the Euler angles, each row must contain three Euler angles "
      "which correspond to each grain orientation. ");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

Eta2EulerAnglesMaterial::Eta2EulerAnglesMaterial(const InputParameters & parameters)
  : Material(parameters),
  _Euler_angles_file_name(getParam<FileName>("Euler_angles_file_name")),
  _op_num(coupledComponents("v")),
  _vals(coupledValues("v")),
  _Euler_angles(_op_num),
  _Euler_angles_mat_prop(declareProperty<RealVectorValue>("Euler_angles"))
{
  for (const auto i : make_range(_op_num)) {

    _Euler_angles[i].zero(); // initialize to zero  

  }

  if (this->isParamValid("Euler_angles_file_name")) {

    assignEulerAngles();
    
  } // otherwise Euler angles are zero
}

void
Eta2EulerAnglesMaterial::computeQpProperties()
{
  RealVectorValue value(0.0, 0.0, 0.0);
  
  for (const auto j : make_range(Moose::dim)) {
    for (unsigned int i = 0; i < _op_num; ++i) {
      value(j) += _Euler_angles[i](j) * ((*_vals[i])[_qp]);
    }
    _Euler_angles_mat_prop[_qp](j) = value(j);
  }
}

void
Eta2EulerAnglesMaterial::assignEulerAngles()
{
  // read in the slip system data from auxiliary text file
  MooseUtils::DelimitedFileReader _reader(_Euler_angles_file_name);
  _reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::ROWS);
  _reader.read();
  
  // add check for the size of the input
  if (_reader.getData().size() != _op_num)
    paramError(
        "_op_num",
        "The number of rows in the Euler angles file should match the number of phase fields.");

  for (const auto i : make_range(_op_num))
  { 
    for (const auto j : index_range(_reader.getData(i))) 
    {
      _Euler_angles[i](j) = _reader.getData(i)[j];
    }
  }
}
