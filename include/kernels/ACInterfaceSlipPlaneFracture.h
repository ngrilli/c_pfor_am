// Nicolò Grilli
// University of Bristol
// 15 Febbraio 2022

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
  
  /// Plane normal to the weak cleavage plane: M in (Clayton & Knap, 2015)
  RealVectorValue _cleavage_plane_normal;

  /// penalty for damage on planes not normal to the weak (favoured) cleavage
  /// plane (Clayton & Knap, 2015)
  const Real _beta_penalty;
  
  // Slip system index to determine slip normal
  // used for the cleavage plane
  const unsigned int _slip_sys_index;
  
  // Slip plane normal to the weak cleavage plane
  // This is assigned to _cleavage_plane_normal
  // based on the slip system index chosen
  const MaterialProperty<std::vector<Real>> & _slip_plane_normals;
  
};
