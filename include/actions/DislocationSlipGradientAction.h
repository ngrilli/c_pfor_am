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

  /// base name for the auxiliary variables
  const std::string _base_name;

  /// Maximum number of active slip systems for the crystalline material being modeled
  const unsigned int _number_slip_systems;
};
