// Nicolo Grilli
// ADSC Singapore
// 13 Ottobre 2020

#pragma once

#include "ComputeElasticityTensorCPGrain.h"

/**
 * ComputeElasticityTensorMelting defines an elasticity tensor material object for crystal plasticity.
 * It is based on the user object GrainPropertyReadFile
 * that allows to assign grains from the euler angles input file
 * to the physical volumes in GMSH
 * Melting is modelled as a degradation of the stiffness tensor
 * when temperature grows above the crystallization temperature
 */
class ComputeElasticityTensorMelting : public ComputeElasticityTensorCPGrain
{
public:
  static InputParameters validParams();

  ComputeElasticityTensorMelting(const InputParameters & parameters);

protected:

  virtual void computeQpElasticityTensor() override;
  
  virtual void melting();
  
  const Real _melting_temperature;
  const Real _crystallisation_temperature;

  /// Stiffness tensor modified by melting
  /// to model laser scanning
  RankFourTensor _Melt_Cijkl;
  
};
