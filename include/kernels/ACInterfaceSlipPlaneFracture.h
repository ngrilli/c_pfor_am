// Nicolò Grilli
// University of Bristol
// 13 Febbraio 2022

// Considers cleavage plane anisotropy in the crack propagation
// the cleavage plane is based on a pre-selected slip plane 

#pragma once

#include "ACInterface.h"

class ACInterfaceSlipPlaneFracture : public ACInterface
{
public:
  static InputParameters validParams();

  ACInterfaceSlipPlaneFracture(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  /// term with beta penalty
  Real betaNablaPsi();

  /// penalty for damage on planes not normal to the weak (favoured) cleavage
  /// plane (Clayton & Knap, 2015)
  const Real _beta_penalty;
  
  // Dislocation velocity value (signed) on all slip systems
  const MaterialProperty<std::vector<Real>> & _dislo_velocity;

  /// Plane normal to the weak cleavage plane: M in (Clayton & Knap, 2015)
  const RealVectorValue _cleavage_plane_normal;
};
