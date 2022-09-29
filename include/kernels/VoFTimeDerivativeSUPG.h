// Nicolò Grilli
// University of Bristol
// 29 Settembre 2022

#pragma once

#include "ADTimeKernelGrad.h"

/**
 * Applies SUPG stabilization to the time derivative.
 * 
 * The difference with LevelSetTimeDerivativeSUPG is the usage of three coupled
 * variables for the components of the velocity.
 */
class VoFTimeDerivativeSUPG : public ADTimeKernelGrad
{
public:
  static InputParameters validParams();

  VoFTimeDerivativeSUPG(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  /// Velocity components
  const VariableValue & _u;
  const VariableValue & _v;
  
  const bool _has_w;
  const VariableValue & _w;
};
