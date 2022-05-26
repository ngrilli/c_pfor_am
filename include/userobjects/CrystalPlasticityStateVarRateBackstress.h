// Nicolò Grilli
// University of Bristol
// 10 Marzo 2022

#pragma once

#include "CrystalPlasticityStateVarRateComponent.h"

/**
 * Evolution rate of the backstress based on Armstrong-Frederick.
 */
class CrystalPlasticityStateVarRateBackstress : public CrystalPlasticityStateVarRateComponent
{
public:
  static InputParameters validParams();

  CrystalPlasticityStateVarRateBackstress(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp,
                                                       std::vector<Real> & val) const;

protected:

  const MaterialProperty<std::vector<Real>> & _mat_prop_slip_rate;

  const MaterialProperty<std::vector<Real>> & _mat_prop_backstress;

  // The backstress parameters read from .i file
  std::vector<Real> _bprops;

};
