// Nicolò Grilli
// University of Bristol
// 19 Gennaio 2023

#pragma once

#include "HEVPFlowRatePowerLawJ2.h"

/**
 * This user object classs
 * Computes flow rate based on power law and
 * Direction based on J2
 * It is similar to HEVPFlowRatePowerLawJ2 but
 * plastic and not visco-plastic, therefore it gives zero
 * plastic strain rate below the yield point
 */
class HEPFlowRatePowerLawJ2 : public HEVPFlowRatePowerLawJ2
{
public:
  static InputParameters validParams();

  HEPFlowRatePowerLawJ2(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;
  virtual bool computeTensorDerivative(unsigned int, const std::string &, RankTwoTensor &) const;

protected:
  Real _threshold_stress_ratio;

};
