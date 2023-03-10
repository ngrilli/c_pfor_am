// Kyprianos Kythreotis
// Nicolò Grilli
// University of Bristol
// 18 Febbraio 2023

#pragma once

#include "ComputeLinearElasticPFFractureStress.h"
#include "MooseEnum.h"

/**
 * Phase-field fracture
 * This class computes the stress and energy contribution for the
 * small strain Linear Elastic formulation of phase field fracture
 * Cyclic model
 */
class ComputeLinearElasticPFFractureCyclic : public ComputeLinearElasticPFFractureStress
{
public:
  static InputParameters validParams();

  ComputeLinearElasticPFFractureCyclic(const InputParameters & parameters);

protected:

  /**
   * Method to split elastic energy based on strain spectral decomposition
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   * Cyclic model
   */
  void computeStrainSpectral(Real & F_pos, Real & F_neg);

  virtual void computeQpStress() override;

  // number of cycles per unit time
  const Real _cycles_per_unit_time;
  
  // Characteristic time for stress decrease in stress history
  const Real _tau_cyclic_stress_history;
  
  // Minimum residual fatigue degradation
  const Real _residual_fatigue_degradation;

  // alpha variable in the cyclic model
  // and its value at the previous time step
  // Miner's rule
  MaterialProperty<Real> & _alpha_cyclic;
  const MaterialProperty<Real> & _alpha_cyclic_old;

  // fatigue degradation function in the cyclic model
  const MaterialProperty<Real> & _fatigue_degradation_old;
  
  // tracks the maximum positive stress eigenvalue
  // to determine the point of the SN curve
  MaterialProperty<Real> & _cyclic_stress_history;
  const MaterialProperty<Real> & _cyclic_stress_history_old;
  
  // Stress level as a function of the number of cycles to failure
  // This is defined using a ParsedMaterial in the input file 
  const MaterialProperty<Real> & _NS_curve_old;

};
