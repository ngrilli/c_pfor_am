// Nicolò Grilli
// University of Bristol
// 29 Settembre 2022

#pragma once

#include "ADKernelValue.h"

/**
 * Advection Kernel for the levelset equation.
 *
 * \psi_i \vec{v} \nabla u,
 * where \vec{v} is the interface velocity that is a set of
 * coupled variables.
 * 
 * The difference with LevelSetAdvection is the usage of three coupled
 * variables for the components of the velocity.
 */
class VolumeOfFluidAdvection : public ADKernelValue
{
public:
  static InputParameters validParams();

  VolumeOfFluidAdvection(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  /// Velocity components
  const VariableValue & _u;
  const VariableValue & _v;
  
  const bool _has_w;
  const VariableValue & _w;
};
