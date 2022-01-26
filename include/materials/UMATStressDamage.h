// Nicolò Grilli
// University of Bristol
// 26 Gennaio 2022

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
  
  // Free energy components and their derivatives
  // calculated by the UMAT
  // are assigned to the corresponding MaterialProperty
  virtual void assignFreeEnergy(Real & F_pos, Real & F_neg);
  
  virtual void computeHistoryVariable(Real & F_pos, Real & F_neg);
  
  /// Variable defining the phase field damage parameter
  const VariableValue & _c;
  
  /// Use current value of history variable
  bool _use_current_hist;
  
  /// History variable that prevents crack healing, declared in this material
  MaterialProperty<Real> & _H;

  /// Old value of history variable
  const MaterialProperty<Real> & _H_old;
  
  /// Material property for elastic energy
  MaterialProperty<Real> & _E;

  /// Derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _dEdc;

  /// Second-order derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _d2Ed2c;
  
  /// Derivative of stress w.r.t damage variable
  MaterialProperty<RankTwoTensor> & _dstress_dc;
  
  /// Second-order derivative of elastic energy w.r.t damage variable and strain
  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;
  
  /// Material property for energetic degradation function
  /// for instance, (1-c)^2 can be used
  const MaterialProperty<Real> & _D;

  /// Derivative of degradation function w.r.t damage variable
  const MaterialProperty<Real> & _dDdc;

  /// Second-order derivative of degradation w.r.t damage variable
  const MaterialProperty<Real> & _d2Dd2c;
  
  const Real _bulk_modulus_ref; // reference bulk modulus for vol/non-vol decomposition

};
