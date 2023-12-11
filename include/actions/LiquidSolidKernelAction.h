// Nicolò Grilli
// Università di Bristol
// 10 Dicembre 2023

#pragma once

#include "InputParameters.h"
#include "Action.h"

/**
 * Action that sets up Reaction, BodyForce, CoupledTanh, ..., Diffusion kernels,
 * for the zeta variable, which is 0 in the liquid phase and 1 in the solid phase.
 */
class LiquidSolidKernelAction : public Action
{
public:
  static InputParameters validParams();

  LiquidSolidKernelAction(const InputParameters & params);

  virtual void act();

protected:

};
