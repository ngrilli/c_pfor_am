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
  /// path of the heat source, x, y, z components
  const Function & _function_x;
  const Function & _function_y;
  const Function & _function_z;

  ADMaterialProperty<Real> & _volumetric_heat;
};
