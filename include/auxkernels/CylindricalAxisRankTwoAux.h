// Nicolò Grilli
// University of Bristol
// 7 Marzo 2023

#pragma once

#include "AuxKernel.h"
#include "RankTwoTensor.h"

/**
 * CylindricalAxisRankTwoAux is designed to take the data in the CylindricalRankTwoTensor material
 * property, for example stress or strain, and output the value for the
 * supplied indices in cylindrical coordinates, where the cylindrical plane axis is
 * along the _axis_vector vector and the center point is defined by _center_point.
 * This is similar to CylindricalRankTwoAux but with an arbitrary axis
 * for the calculation of the cylindrical coordinates.
 */

class CylindricalAxisRankTwoAux : public AuxKernel
{
public:
  static InputParameters validParams();

  CylindricalAxisRankTwoAux(const InputParameters & parameters);
  virtual ~CylindricalAxisRankTwoAux() {}

protected:
  virtual Real computeValue();
  const MaterialProperty<RankTwoTensor> & _tensor;
  const unsigned int _i;
  const unsigned int _j;
  const Point _center_point;
  
  /// User-input value of the vector parallel to the axis with respect to which
  /// cylindrical coordinates are calculated
  RealVectorValue _axis_vector;
};
