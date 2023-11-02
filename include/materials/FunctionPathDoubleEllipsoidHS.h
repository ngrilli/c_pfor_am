// Nicolò Grilli
// Valentin Godard
// University of Bristol
// 2 Novembre 2023

#pragma once

#include "Material.h"

class Function;

/**
 * Double ellipsoid heat source distribution based on:
 * O. Mokrov, M. Simon, A. Schiebahn, U. Reisgen,
 * A fine modification of the double ellipsoid heat source,
 * Mathematical Modelling of Weld Phenomena 12
 */
class FunctionPathDoubleEllipsoidHS : public Material
{
public:
  static InputParameters validParams();

  FunctionPathDoubleEllipsoidHS(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// power
  const Real _P;
  /// process efficienty
  const Real _eta;
  
  /// ellipsoid axes along the scanning direction: front and back
  const Real _a_f;
  const Real _a_r;
  
  /// ellipsoid axes perpendicular to the scanning direction
  const Real _b;
  const Real _c;
  
  /// scaling factor for the front part
  const Real _f_f;
  /// scaling factor for the rear part
  const Real _f_r;
  
  /// path of the heat source, x, y, z components
  const Function & _function_x;
  const Function & _function_y;
  const Function & _function_z;

  ADMaterialProperty<Real> & _volumetric_heat;
};
