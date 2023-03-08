// Nicolò Grilli
// University of Bristol
// 7 Marzo 2023

#include "CylindricalAxisRankTwoAux.h"

registerMooseObject("TensorMechanicsApp", CylindricalAxisRankTwoAux);

InputParameters
CylindricalAxisRankTwoAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Takes RankTwoTensor material and outputs component in cylindrical coordinates. "
      "This is similar to CylindricalRankTwoAux but with an arbitrary axis "
      "for the calculation of the cylindrical coordinates. ");
  params.addRequiredParam<MaterialPropertyName>("rank_two_tensor",
                                                "The rank two material tensor name");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i",
      "index_i >= 0 & index_i <= 2",
      "The index i of ij for the tensor to output (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j",
      "index_j >= 0 & index_j <= 2",
      "The index j of ij for the tensor to output (0, 1, 2)");
  params.addRequiredParam<Point>("center_point",
                                 "Location of the center point of the cylindrical coordinates");
  params.addRequiredParam<RealVectorValue>("axis_vector", "The vector parallel to the axis with respect to which "
                                                          "cylindrical coordinates are calculated. ");
  return params;
}

CylindricalAxisRankTwoAux::CylindricalAxisRankTwoAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _tensor(getMaterialProperty<RankTwoTensor>("rank_two_tensor")),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j")),
    _center_point(getParam<Point>("center_point")),
    _axis_vector(getParam<RealVectorValue>("axis_vector"))
{
  // Normalize axis vector
  if (_axis_vector.norm() == 0)
    mooseError("CylindricalAxisRankTwoAux: _axis_vector must not have zero length");
  else
    _axis_vector /= _axis_vector.norm();  
}

Real
CylindricalAxisRankTwoAux::computeValue()
{		
  // Distance between current point and _center_point
  Point loc_from_center = _q_point[_qp] - _center_point;
  
  // Scalar projection of loc_from_center along _axis_vector
  Real projection_along_axis_vector = 0.0;
  
  // Projected distance vector from the present coordinates to _center_point
  RealVectorValue projected_loc_from_center;
  
  RealVectorValue x_axis(1.0, 0.0, 0.0);
  RealVectorValue y_axis(0.0, 1.0, 0.0);
  RealVectorValue z_axis(0.0, 0.0, 1.0);
  
  // Axes defining a plane perpendicular to _axis_vector
  RealVectorValue xprime_axis, yprime_axis;
  
  // xprime and yprime coordinates of projected_loc_from_center
  Real xprime_loc_from_center = 0.0;
  Real yprime_loc_from_center = 0.0;
  
  // Projection of the z axis onto _axis_vector
  Real z_projection_on_axis_vector = 0.0;
  
  // Input tensor in the primed reference system
  // and rotation matrix for the transformation
  RankTwoTensor primed_tensor;
  RankTwoTensor Primed_to_Sample;
  
  for (const auto j : make_range(LIBMESH_DIM)) {
	  
    projection_along_axis_vector += loc_from_center(j) * _axis_vector(j);
  
  }
  
  // Project loc_from_center on a plane perpendicular to _axis_vector
  for (const auto j : make_range(LIBMESH_DIM)) {
	  
    projected_loc_from_center(j) = loc_from_center(j) - projection_along_axis_vector * _axis_vector(j);
	  
  } 
  
  // Rotate the x and y axes to obtain new axes xprime and yprime 
  // to calculate the angle
  // theta for cylindrical coordinates conversion
  
  z_projection_on_axis_vector = _axis_vector(2);
  
  if (z_projection_on_axis_vector < 0.999) {
	  
    xprime_axis = z_axis.cross(_axis_vector);
    xprime_axis /= xprime_axis.norm();
    yprime_axis = _axis_vector.cross(xprime_axis);
	  
  } else { // Case in which _axis_vector is almost like z axis

    xprime_axis = x_axis;
    yprime_axis = y_axis;

  }

  // Calculate cylindrical coordinates in the new primed reference system
  
  for (const auto j : make_range(LIBMESH_DIM)) {
	  
    xprime_loc_from_center += projected_loc_from_center(j) * xprime_axis(j);
    yprime_loc_from_center += projected_loc_from_center(j) * yprime_axis(j);  
	  
  }
  
  // Build a rotation matrix to transform the coordinates of a vector in the primed system
  // into the coordinates on the sample (x,y,z) system
  
  for (const auto j : make_range(LIBMESH_DIM)) {
	  
    Primed_to_Sample(j,0) = xprime_axis(j);
    Primed_to_Sample(j,1) = yprime_axis(j);	  
    Primed_to_Sample(j,2) = _axis_vector(j); 	  
	  
  }

  // Transform the tensor into the primed system
  
  primed_tensor = Primed_to_Sample.transpose() * _tensor[_qp] * Primed_to_Sample;

  // Transform into cylindrical coordinates

  Real theta = std::atan2(yprime_loc_from_center, xprime_loc_from_center);
  RankTwoTensor R;
  R(0, 0) = std::cos(theta);
  R(0, 1) = std::sin(theta);
  R(1, 0) = -std::sin(theta);
  R(1, 1) = std::cos(theta);
  R(2, 2) = 1.0;

  RankTwoTensor rotated_tensor = R * primed_tensor * R.transpose();

  return rotated_tensor(_i, _j);
}
