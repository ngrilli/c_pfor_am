// Nicolò Grilli
// Università di Bristol
// 22 Novembre 2025

// This is a modification of FVAdvection
// where the velocity vector is imported from a material class
// that calculates dislocation velocity 
// based on crystallographic directions and stress

#pragma once

#include "FVFluxKernel.h"

class DislocationFVAdvection : public FVFluxKernel
{
public:
  static InputParameters validParams();
  DislocationFVAdvection(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _velocity;

  /// The interpolation method to use for the advected quantity
  Moose::FV::InterpMethod _advected_interp_method;
  
  /// Edge slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _edge_slip_direction;

  /// Screw slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
  /// Slip system index to determine slip direction
  const unsigned int _slip_sys_index;
};
