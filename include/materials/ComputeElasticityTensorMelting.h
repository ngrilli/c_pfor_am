// Nicolo Grilli
// National University of Singapore
// 10 Novembre 2020

#pragma once

#include "ComputeElasticityTensorCPGrain.h"

/**
 * ComputeElasticityTensorMelting defines an elasticity tensor material object for crystal plasticity.
 * It is based on the user object GrainPropertyReadFile
 * that allows to assign grains from the euler angles input file
 * to the physical volumes in GMSH
 * Melting is considered: stiffness is degraded when the
 * temperature increases above melting or below gas temperature.
 */
class ComputeElasticityTensorMelting : public ComputeElasticityTensorCPGrain
{
public:
  static InputParameters validParams();

  ComputeElasticityTensorMelting(const InputParameters & parameters);

protected:

  virtual void computeQpElasticityTensor() override;
  
  virtual void melting();
  
  const Real _melting_temperature_high;
  const Real _melting_temperature_low;
  const Real _gas_temperature_high;
  const Real _gas_temperature_low;
  
  /// Residual stiffness of gas and molten pool (percent)
  const Real _residual_stiffness;

  /// Stiffness tensor modified by melting
  /// to model laser scanning
  RankFourTensor _Melt_Cijkl;
  
};
