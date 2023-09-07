// Nicolò Grilli
// Università di Bristol
// 3 Settembre 2023

#pragma once

#include "Material.h"

/**
 * This class increases a material property from 0 to 1 if the element is inside an ellipsoid
 */
class FunctionPathEllipsoidMaterial : public Material
{
public:
  static InputParameters validParams();

  FunctionPathEllipsoidMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeQpProperties() override;

  /// Base name prepended to material property name
  const std::string _base_name;
  
  // The default is 0 to 1, but in general the level set will be transformed from
  // _low_level_set_var to _high_level_set_var
  const Real _low_level_set_var;
  const Real _high_level_set_var;
  
  /// transverse ellipsoid axe
  const Real _rx;
  /// depth ellipsoid axe
  const Real _ry;
  /// longitudinal ellipsoid axe
  const Real _rz;

  /// path of the heat source, x, y, z components
  const Function & _function_x;
  const Function & _function_y;
  const Function & _function_z;
  
  // Threshold value of the ellipsoid function
  // that activates the level set.
  const Real _level_set_activation_threshold;

  MaterialProperty<Real> & _level_set;
  const MaterialProperty<Real> & _level_set_old;
  
};
