// Nicolò Grilli
// University of Bristol
// 17 Marzo 2022

#pragma once

#include "CrystalPlasticityStateVariable.h"

/**
 * Crystal plasticity state variable userobject class.
 * It can handle state variables
 * that can become positive or negative during time.
 */
class CrystalPlasticityPosNegStateVariable : public CrystalPlasticityStateVariable
{
public:
  static InputParameters validParams();

  CrystalPlasticityPosNegStateVariable(const InputParameters & parameters);

  virtual bool updateStateVariable(unsigned int qp,
                                   Real dt,
                                   std::vector<Real> & val,
                                   std::vector<Real> & val_old) const;

};
