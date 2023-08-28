// Nicolò Grilli
// Università di Bristol
// 28 Agosto 2023

#pragma once

#include "InputParameters.h"
#include "PolycrystalVariablesAction.h"

/**
 * Automatically generates all variables to model a polycrystal with op_num orderparameters
 * and zeta solidification variable
 */
class SolidificationVariablesAction : public PolycrystalVariablesAction
{
public:
  static InputParameters validParams();

  SolidificationVariablesAction(const InputParameters & params);

  virtual void act();

protected:
  // These variables need to be redefined because they are private
  // in the parent class
  const unsigned int _op_num_solidification;
  const std::string _var_name_base_solidification;

};
