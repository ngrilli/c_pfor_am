// Nicolò Grilli
// Università di Bristol
// 26 Agosto 2023

#pragma once

#include "PolycrystalKernelAction.h"

/**
 * Action that sets up ACGrGrPoly, GrainSolidification, ACInterface, TimeDerivative, and ACGBPoly
 * kernels.
 */
class PolycrystalSolidificationKernelAction : public PolycrystalKernelAction
{
public:
  static InputParameters validParams();

  PolycrystalSolidificationKernelAction(const InputParameters & params);

  virtual void act();

protected:

};
