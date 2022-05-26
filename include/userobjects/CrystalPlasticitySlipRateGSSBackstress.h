// Nicolò Grilli
// University of Bristol
// 17 Marzo 2022

#pragma once

#include "CrystalPlasticitySlipRateGSS.h"
#include "RankTwoTensor.h"

/**
 * Phenomenological constitutive model slip rate userobject class.
 * Backstress included
 */
class CrystalPlasticitySlipRateGSSBackstress : public CrystalPlasticitySlipRateGSS
{
public:
  static InputParameters validParams();

  CrystalPlasticitySlipRateGSSBackstress(const InputParameters & parameters);

  virtual bool calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual bool calcSlipRateDerivative(unsigned int qp, Real /*dt*/, std::vector<Real> & val) const;

protected:

  // Backstress read from state variable user object 
  const MaterialProperty<std::vector<Real>> & _mat_prop_backstress;

};
