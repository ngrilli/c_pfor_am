// Nicolò Grilli
// University of Bristol
// 13 Febbraio 2022

#pragma once

#include "CrystalPlasticitySlipRateGSS.h"
#include "RankTwoTensor.h"

/**
 * Phenomenological constitutive model slip rate userobject class.
 * this is the same as CrystalPlasticitySlipRateGSS
 * but can output the slip plane normal using the
 * function calcFlowDirection
 */
class CrystalPlasticitySlipRateCleavage : public CrystalPlasticitySlipRateGSS
{
public:
  static InputParameters validParams();

  CrystalPlasticitySlipRateCleavage(const InputParameters & parameters);

  virtual void calcFlowDirection(unsigned int qp,
                                 std::vector<RankTwoTensor> & flow_direction,
								 std::vector<Real> & slip_plane_normals) const;
  
};
