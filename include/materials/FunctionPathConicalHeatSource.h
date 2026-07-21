// Nicolò Grilli
// Zhuohao Song
// Università di Bristol

#pragma once

#include "Material.h"

class Function;

/**
 * Sum between two ellipsoidal heat source distribution 
 * and a conical heat source distribution.
 * J. Wang, N.O. Larrosa, C. Jacquemoud, C.E. Truman
 * Forward and inverse surrogate modelling of welding residual stress in electron beam welded 316L steel using multilayer perceptron networks
 * International Journal of Pressure Vessels and Piping, Volume 222, Part 2 (2026) 105826
 * https://www.sciencedirect.com/science/article/pii/S0308016126000827
 */
class FunctionPathConicalHeatSource : public Material
{
public:
  static InputParameters validParams();

  FunctionPathConicalHeatSource(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// power
  const Real _P_ellipsoid_1;
  const Real _P_ellipsoid_2;
  const Real _P_conical;
  /// process efficienty
  const Real _eta_ellipsoid_1;
  const Real _eta_ellipsoid_2;
  const Real _eta_conical;
  /// ellipsoidal axes
  const Real _rx_ellipsoid_1;
  const Real _ry_ellipsoid_1;
  const Real _rz_ellipsoid_1;
  const Real _rx_ellipsoid_2;
  const Real _ry_ellipsoid_2;
  const Real _rz_ellipsoid_2;
  /// conical axes
  const Real _rx_conical;
  const Real _ry_conical;
  const Real _rz_conical;
  const Real _cone_radius_decrease;
  /// scaling factor
  const Real _f_ellipsoid_1;
  const Real _f_ellipsoid_2;
  const Real _f_conical;
  /// path of the heat sources, x, y, z components
  const Function & _function_x_ellipsoid_1;
  const Function & _function_y_ellipsoid_1;
  const Function & _function_z_ellipsoid_1;
  const Function & _function_x_ellipsoid_2;
  const Function & _function_y_ellipsoid_2;
  const Function & _function_z_ellipsoid_2;
  const Function & _function_x_conical;
  const Function & _function_y_conical;
  const Function & _function_z_conical;

  ADMaterialProperty<Real> & _volumetric_heat;
};
