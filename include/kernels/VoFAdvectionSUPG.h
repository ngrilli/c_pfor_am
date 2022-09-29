// Nicolò Grilli
// University of Bristol
// 29 Settembre 2022

#pragma once

#include "ADKernelGrad.h"

/**
 * SUPG stabilization for the advection portion of the level set equation.
 * 
 * The difference with LevelSetAdvectionSUPG is the usage of three coupled
 * variables for the components of the velocity.
 */
class VoFAdvectionSUPG : public ADKernelGrad
{
public:
  static InputParameters validParams();

  VoFAdvectionSUPG(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  /// Velocity components
  const VariableValue & _u;
  const VariableValue & _v;
  
  const bool _has_w;
  const VariableValue & _w;
};
