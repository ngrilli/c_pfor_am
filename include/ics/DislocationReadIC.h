// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// Univerisity of Bristol
// 25 Agosto 2022

#pragma once

#include "InitialCondition.h"
#include "PropertyReadFile.h"

/**
 * Initialize dislocation loops read from file
 * File must contain four columns with:
 * rho_tot rho_e rho_s q
 */
class DislocationReadIC : public InitialCondition
{
public:
  static InputParameters validParams();

  DislocationReadIC(const InputParameters & parameters);

protected:

  virtual Real value(const Point & p) override;

  /**
   * Element property read user object
   * Presently used to read the four dislocation density variables 
   */
  const PropertyReadFile * const _read_dislocation_user_object;

  // Type of variable for the four dislocation density variables
  // rhotot rhoedgegnd rhoscrewgnd qtot
  MooseEnum _variable_type;

};
