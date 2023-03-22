// Fernando Valiente Dies
// ANSTO
// Nicolò Grilli
// University of Bristol
// 22 Marzo 2023

#pragma once

#include "AuxKernel.h"

class Function;

/**
 * This AuxKernel increases an AuxVariable from 0 to 1 if the qp is inside
 * an ellipsoid that is moving according to paths defined by velocity and initial coordinates.
 * It can be applied to a level set variable
 * to simulate the material deposition during wire arc additive manufacturing (WAAM)
 * together with ActDeactElementsCoupled.
 * It uses a postprocessor to start the next path to simulate a temperature controlled WAAM deposition process
 */
class VelocityEllipsoidAux : public AuxKernel
{
public:
  static InputParameters validParams();

  VelocityEllipsoidAux(const InputParameters & parameters);

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

  /// Velocity vector of the material deposition source
  const RealVectorValue _velocity;

  /// Initial values of the coordinates of the material deposition source
  /// Every time the postprocessor condition is satisfied, 
  /// the heat source is moved to the next set of coordinates
  const std::vector<Real> _init_x_coords;
  const std::vector<Real> _init_y_coords;
  const std::vector<Real> _init_z_coords;
  
  /// Postprocess with temperature value
  /// it provides the condition based on which the material deposition source
  /// is moved to the next set of initial coordinates
  const PostprocessorValue & _temperature_pp;
  
  
  
  
  
  
  // Threshold value of the ellipsoid function
  // that activates the level set.
  const Real _level_set_activation_threshold;

};
