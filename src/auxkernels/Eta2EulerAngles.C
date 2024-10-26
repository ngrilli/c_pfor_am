// Nicolò Grilli
// Università di Bristol
// Fernando Valiente Dies
// ANSTO
// 26 Ottobre 2024

#include "Eta2EulerAngles.h"

registerMooseObject("c_pfor_amApp", Eta2EulerAngles);

InputParameters
Eta2EulerAngles::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Output euler angles from eta_i phase fields");
  MooseEnum euler_angles("phi1 Phi phi2");
  params.addRequiredParam<MooseEnum>("output_euler_angle", euler_angles, "Euler angle to output");
  params.addParam<FileName>(
      "Euler_angles_file_name","",
      "Name of the file containing the Euler angles, each row must contain three Euler angles "
      "which correspond to each grain orientation. ");
  params.addRequiredCoupledVar("v",
                               "Array of coupled order parameter names for the eta variables");
  return params;
}

Eta2EulerAngles::Eta2EulerAngles(const InputParameters & parameters)
  : AuxKernel(parameters),
    _output_euler_angle(getParam<MooseEnum>("output_euler_angle")),
    _Euler_angles_file_name(getParam<FileName>("Euler_angles_file_name")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _Euler_angles(_op_num)
{
  for (const auto i : make_range(_op_num)) {

    _Euler_angles[i].zero(); // initialize to zero  

  }

  if (this->isParamValid("Euler_angles_file_name")) {

    assignEulerAngles();
    
  } // otherwise Euler angles are zero
}

void
Eta2EulerAngles::assignEulerAngles()
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

void
Eta2EulerAngles::precalculateValue()
{
  _value = 0.0;
  
  for (unsigned int i = 0; i < _op_num; ++i)
    _value += _Euler_angles[i](_output_euler_angle) * (*_vals[i])[_qp];
}

Real
Eta2EulerAngles::computeValue()
{
  return _value;
}
