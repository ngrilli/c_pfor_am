// Nicolò Grilli
// Università di Bristol
// 23 Dicembre 2024

#pragma once

#include "Action.h"

/**
 * Action that sets up slip rate auxkernels for the dislocation slip gradient model
 */

class DislocationSlipGradientAction : public Action
{
public:
  static InputParameters validParams();

  DislocationSlipGradientAction(const InputParameters & params);

  virtual void act();

protected:

  /// Maximum number of active slip systems for the crystalline material being modeled
  const unsigned int _number_slip_systems;
};
