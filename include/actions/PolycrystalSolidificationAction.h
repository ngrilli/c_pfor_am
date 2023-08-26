// Nicolò Grilli
// Università di Bristol
// 26 Agosto 2023

#pragma once

#include "PolycrystalKernelAction.h"

/**
 * Action that sets up ACGrGrPoly, GrainSolidification, ACInterface, TimeDerivative, and ACGBPoly
 * kernels.
 */
class PolycrystalSolidificationAction : public PolycrystalKernelAction
{
public:
  static InputParameters validParams();

  PolycrystalSolidificationAction(const InputParameters & params);

  virtual void act();

protected:

};
