// Nicolò Grilli
// Università di Bristol
// 28 Febbraio 2024

#pragma once

#include "Action.h"

/**
 * Action that sets up slip rate auxkernels to visualize dislocation densities
 */

class DislocationVariablesAction : public Action
{
public:
  static InputParameters validParams();

  DislocationVariablesAction(const InputParameters & params);

  virtual void act();

protected:

  /// base name for the auxiliary variables
  const std::string _base_name;

  /// Maximum number of active slip systems for the crystalline material being modeled
  const unsigned int _number_slip_systems;
  
  /// Flags to select which dislocation densities to visualize
  const bool _rho_ssd;
  const bool _rho_gnd_edge;
  const bool _rho_gnd_screw;
  const bool _slip_resistance;
};
