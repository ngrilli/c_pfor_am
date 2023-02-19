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

  // alpha variable in the cyclic model
  // and its value at the previous time step
  MaterialProperty<Real> & _alpha_cyclic;
  const MaterialProperty<Real> & _alpha_cyclic_old;

  // fatigue degradation function in the cyclic model
  MaterialProperty<Real> & _fatigue_degradation;


};
