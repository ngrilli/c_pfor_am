// Nicolò Grilli
// University of Bristol
// 26 Settembre 2021

#pragma once

#include "FiniteStrainUObasedCP.h"

// Are these needed?
#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "CrystalPlasticityStateVariable.h"
#include "CrystalPlasticityStateVarRateComponent.h"

/**
 * FiniteStrainUObasedCPDamage is the coupling between
 * crystal plasticity and damage
 */
class FiniteStrainUObasedCPDamage : public FiniteStrainUObasedCP
{
public:
  static InputParameters validParams();

  FiniteStrainUObasedCPDamage(const InputParameters & parameters);

protected:
  /**
   * calculate stress residual.
   * Damage is added to the stress calculation
   */
  virtual void calcResidual();
  
  /**
   * Method to split elastic energy based on strain volumetric/deviatoric decomposition
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   * @param ee elastic Green-Lagrange strain
   */
  virtual void computeStrainSpectral(Real & F_pos, Real & F_neg, 
                                     RankTwoTensor & ee, RankTwoTensor & pk2_new);

  // compute history variable and assign to _E
  // which is used by the fracture model for damage growth
  // Damage grows only because of the positive part of the elastic energy F_pos									 
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

};
