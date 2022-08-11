// Nicolò Grilli
// University of Bristol
// 7 Agosto 2022

#pragma once

#include "AuxKernel.h"

class Function;

/**
 * This AuxKernel increases an AuxVariable from 0 to 1 if the qp is inside
 * an ellipsoid that is moving according to a path provided through a function.
 * It can be applied to a level set variable
 * to simulate the material deposition during wire arc additive manufacturing
 * together with ActDeactElementsCoupled.
 */
class FunctionPathEllipsoidAux : public AuxKernel
{
public:
  static InputParameters validParams();

  FunctionPathEllipsoidAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  
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

};
