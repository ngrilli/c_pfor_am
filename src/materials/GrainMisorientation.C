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
    _misorientation(declareProperty<Real>("misorientation"))
{
}

void
GrainMisorientation::initQpStatefulProperties()
{

}

void
GrainMisorientation::computeQpProperties()
{
  //computeQpStress();


}
