// Nicolò Grilli
// Università di Bristol
// 26 Agosto 2023

#pragma once

#include "PolycrystalKernelAction.h"

/**
 * Action that sets up ACGrGrPoly, GrainSolidification, ACInterface, TimeDerivative, and ACGBPoly
 * kernels, based on (5) in:
 * Min Yang, Lu Wang and Wentao Yan
 * Phase-ﬁeld modeling of grain evolutions in additive
 * manufacturing from nucleation, growth, to coarsening
 * npj Computational Materials volume 7, 56 (2021)
 */
class PolycrystalSolidificationKernelAction : public PolycrystalKernelAction
{
public:
  static InputParameters validParams();

  PolycrystalSolidificationKernelAction(const InputParameters & params);

  virtual void act();

protected:

};
