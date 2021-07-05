// Nicolò Grilli
// University of Bristol
// 5 Luglio 2021

#pragma once

#include "FVCoupledFluxKernel.h"

class FVCoupledAdvection : public FVCoupledFluxKernel
{
public:
  static InputParameters validParams();
  FVCoupledAdvection(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

  const RealVectorValue _velocity;

  /// The interpolation method to use for the advected quantity
  Moose::FV::InterpMethod _advected_interp_method;
};
