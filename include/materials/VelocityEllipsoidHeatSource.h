// Fernando Valiente Dies & Nicolò Grilli
// University of Sydney & University of Bristol
// 26 Novembre 2022

// Double ellipsoid volumetric source heat, the motion is determined
// by input velocity, starting positions and a postprocessor.

#pragma once

#include "Material.h"

//class Function;  // is this necessary?

/**
 * Double ellipsoid heat source distribution.
 */
class VelocityEllipsoidHeatSource : public Material
{
public:
  static InputParameters validParams();

  VelocityEllipsoidHeatSource(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  
  // Check if the postprocessor temperature condition is satisfied
  // and change the initial coordinates and scan time
  virtual void checkPPcondition();

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
  const RealVectorValue _velocity;
  
  /// Postprocessor temperature at the previous timestep
  Real _previous_pp_temperature;
  
  /// Initial values of the coordinates of the heat source
  /// Every time the postprocessor condition is satisfied, 
  /// the heat source is moved to the next set of coordinates
  const std::vector<Real> _init_x_coords;
  const std::vector<Real> _init_y_coords;
  const std::vector<Real> _init_z_coords;
  
  /// Postprocess with temperature value
  /// it provides the condition based on which the heat source
  /// is moved to the next set of initial coordinates
  const PostprocessorValue & _temperature_pp;
  
  /// variables to store the coordinates 
  /// of the center of the heat source
  Real _x_coord;
  Real _y_coord;
  Real _z_coord;
  
  /// _t_scan tracks the simulation time at which a new
  /// scan begins after the condition based on the postprocessor
  /// changes the coordinates of the heat source
  Real _t_scan;
  
  /// Total length during each scan
  /// After the laser has travelled this length
  /// it is switched off
  const std::vector<Real> _scan_length;
  
  /// When the temperature provided by the postprocessor decreases
  /// below this threshold, the heat source is moved to the next
  /// set of coordinates
  const Real _threshold_temperature;
  
  /// Heat source track index
  int _n_track; 

  ADMaterialProperty<Real> & _volumetric_heat;
};
