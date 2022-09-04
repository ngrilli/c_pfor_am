// Nicolò Grilli
// University of Bristol
// 4 Settembre 2022

#pragma once

#include "Material.h"

/**
 * Calculate misorientation between neighbouring elements
 * of a structured mesh
 */
class GrainMisorientation : public Material
{
public:
  static InputParameters validParams();

  GrainMisorientation(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Misorientation between neighboring grains
  MaterialProperty<Real> & _misorientation;
};
