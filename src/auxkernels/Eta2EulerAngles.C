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
  params.addClassDescription("Output euler angles from eta_i phase fields.");
  MooseEnum euler_angles("phi1 Phi phi2");
  params.addRequiredParam<MooseEnum>("output_euler_angle", euler_angles, "Euler angle to output");
  params.addParam<FileName>(
      "Euler_angles_file_name","",
      "Name of the file containing the Euler angles, each row must contain three Euler angles "
      "which correspond to each grain orientation. ");
  return params;
}

Eta2EulerAngles::Eta2EulerAngles(const InputParameters & parameters)
  : AuxKernel(parameters),
    _output_euler_angle(getParam<MooseEnum>("output_euler_angle")),
    _Euler_angles_file_name(getParam<FileName>("Euler_angles_file_name"))
{






}

void
Eta2EulerAngles::precalculateValue()
{
  // ID of unique grain at current point
//  const auto grain_id =
//      _grain_tracker.getEntityValue((isNodal() ? _current_node->id() : _current_elem->id()),
//                                    FeatureFloodCount::FieldType::UNIQUE_REGION,
//                                    0);

  // Recover euler angles for current grain
//  RealVectorValue angles;
//  if (grain_id >= 0)
//    angles = _euler.getEulerAngles(grain_id);

  // Return specific euler angle
//  _value = angles(_output_euler_angle);
  _value = 0.0;
}

Real
Eta2EulerAngles::computeValue()
{
  return _value;
}
