// Nicolò Grilli
// 28 Ottobre 2020
// National University of Singapore

#pragma once

#include "AuxKernel.h"

// Compute temperature due to laser scan during SLM from CFD simulation

// Forward declarations

class LaserScanAux : public AuxKernel
{
public:
  static InputParameters validParams();

  LaserScanAux(const InputParameters & parameters);
  virtual ~LaserScanAux() {}

protected:
  virtual Real computeValue();

  /// Base name of the material system used to calculate the temperature
  const std::string _base_name;
  
  /// Laser scan velocity magnitude, along positive x coordinates
  const Real _scan_velocity;
  
  std::vector<Real> _laser_init_coord;

};
