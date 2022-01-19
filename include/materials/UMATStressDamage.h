// Nicolò Grilli
// University of Bristol
// 19 Gennaio 2022

#pragma once

#include "AbaqusUMATStress.h"
#include "DynamicLibraryLoader.h"

/**
 * Coupling material to use Abaqus UMAT models in MOOSE
 * with phase field damage
 */
class UMATStressDamage : public AbaqusUMATStress
{
public:
  static InputParameters validParams();

  UMATStressDamage(const InputParameters & parameters);

protected:

  void computeQpStress() override;
  
  /// Variable defining the phase field damage parameter
  const VariableValue & _c;

};
