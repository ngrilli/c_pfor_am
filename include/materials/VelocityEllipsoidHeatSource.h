// Nicolò Grilli
// University of Bristol
// 6 Novembre 2022

// Double ellipsoid volumetric source heat, the motion is determined
// by input velocity, starting positions and a postprocessor.

#pragma once

#include "Material.h"

class Function;  // is this necessary?

/**
 * Double ellipsoid heat source distribution.
 */
class VelocityEllipsoidHeatSource : public Material
{
public:
  static InputParameters validParams();

  VelocityEllipsoidHeatSource(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// power
  const Real _P;
  /// process efficienty
  const Real _eta;
  /// transverse ellipsoid axe
  const Real _rx;
  /// depth ellipsoid axe
  const Real _ry;
  /// longitudinal ellipsoid axe
  const Real _rz;
  /// scaling factor
  const Real _f;
  
  /// Scanning speed vector
  RealVectorValue _velocity;
  
  /// Initial values of the coordinates of the heat source
  /// Every time the postprocessor condition is satisfied, 
  /// the heat source is moved to the next set of coordinates
  std::vector<Real> _init_x_coords;
  std::vector<Real> _init_y_coords;
  std::vector<Real> _init_z_coords;
  
  
  /// path of the heat source, x, y, z components
//  const Function & _function_x;
//  const Function & _function_y;
//  const Function & _function_z;

  ADMaterialProperty<Real> & _volumetric_heat;
};
