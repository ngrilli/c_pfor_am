// Nicolò Grilli
// University of Bristol
// Daijun Hu
// National University of Singapore
// 9 Settembre 2022

#pragma once

#include "ComputeMultipleCrystalPlasticityStress.h"

#include "CrystalPlasticityDislocationUpdateBase.h"
#include "ComputeCrystalPlasticityEigenstrainBase.h"
#include "ElementPropertyReadFile.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"

/**
 * ComputeCrystalPlasticityStressDamage (used together with CrystalPlasticityDislocationUpdateBase)
 * uses the multiplicative decomposition of the deformation gradient and solves the PK2 stress
 * residual equation at the intermediate configuration to evolve the material state. The internal
 * variables are updated using an iterative predictor-corrector algorithm. Backward Euler
 * integration rule is used for the rate equations.
 * Thermal eigenstrain is included
 *
 * The only difference between this class and ComputeMultipleCrystalPlasticityStress
 * is that the _models variable here is an array of CrystalPlasticityDislocationUpdateBase
 * instead of CrystalPlasticityStressUpdateBase.
 *
 * This material model is coupled with phase field damage.
 *
 * Damage formulation is based on:
 * Nicolo Grilli and Marisol Koslowski
 * The effect of crystal anisotropy and plastic response
 * on the dynamic fracture of energetic materials
 * Journal of Applied Physics 126, 155101 (2019).
 * Plastic work is added to the positive part of the free energy that causes damage
 */
class ComputeCrystalPlasticityStressDamage : public ComputeMultipleCrystalPlasticityStress
{
public:
  static InputParameters validParams();

  ComputeCrystalPlasticityStressDamage(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /**
   * Updates the stress (PK2) at a quadrature point by calling constiutive
   * relationship as defined in a child class
   * Solves stress residual equation using Newton - Rhapson: Updates slip
   * system resistances iteratively
   */
  virtual void updateStress(RankTwoTensor & cauchy_stress, RankFourTensor & jacobian_mult) override;

  /**
   * initializes the stateful properties such as PK2 stress, resolved shear
   * stress, plastic deformation gradient, slip system resistances, etc.
   * This class is often overwritten by inherting classes.
   */
  virtual void initQpStatefulProperties() override;

  /**
   * Reset the PK2 stress and the inverse deformation gradient to old values and
   * provide an interface for inheriting classes to reset material properties
   */
  void preSolveQp();

  /**
   * Solve the stress and internal state variables (e.g. slip increment,
   * slip system resistance) at each qp points
   */
  void solveQp();

  /**
   * Save the final stress and internal variable values after the iterative solve.
   */
  void postSolveQp(RankTwoTensor & stress_new, RankFourTensor & jacobian_mult);
  
  /**
   * Plastic work updated according to
   * equation (17) in:
   * Elastic plastic deformation at finite strains
   * E. H. Lee 1968,
   * Stanford University technical report AD678483
   */
  void updatePlasticWork();
  
  /**
   * Solves the internal variables stress as a function of the slip specified
   * by the constitutive model defined in the inheriting class
   */
  void solveStateVariables();
  
  /**
   * solves for stress, updates plastic deformation gradient.
   */
  void solveStress();
  
  /**
   * Calls the residual and jacobian functions used in the stress update
   * algorithm.
   */
  void calculateResidualAndJacobian();

  /**
   * Calculate stress residual as the difference between the stored material
   * property PK2 stress and the elastic PK2 stress calculated from the
   * constitutively defined equivalent_slip_increment.
   * The equivalent_slip_increment is passed in as an input arguement.
   */
  void calculateResidual();

  /**
   * Calculates the jacobian as
   * $\mathbf{J} = \mathbf{I} - \mathbf{C} \frac{d\mathbf{E}^e}{d\mathbf{F}^e}
   * \frac{d\mathbf{F}^e}{d\mathbf{F}^P^{-1}} \frac{d\mathbf{F}^P^{-1}}{d\mathbf{PK2}}$
   */
  void calculateJacobian();
  
  /**
   * Method to split elastic energy based on strain volumetric/non-volumetric decomposition
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   * @param ee elastic Green-Lagrange strain
   * @param ce right Cauchy–Green deformation tensor (elastic)
   * @param pk2_new second Piola-Kirchhoff stress
   */
  virtual void computeStrainVolumetric(Real & F_pos, Real & F_neg, 
                                       RankTwoTensor & ee, RankTwoTensor & ce, 
									   RankTwoTensor & pk2_new);
									   
  /// compute history variable and assign to _E
  /// which is used by the fracture model for damage growth
  /// Damage grows only because of the positive part of the elastic energy F_pos									 
  virtual void computeHistoryVariable(Real & F_pos, Real & F_neg);

  ///@{Calculates the tangent moduli for use as a preconditioner, using the elastic or elastic-plastic option as specified by the user
  void elastoPlasticTangentModuli(RankFourTensor & jacobian_mult);
  ///@}

  /// The user supplied cyrstal plasticity consititutive models
  std::vector<CrystalPlasticityDislocationUpdateBase *> _dislocation_models;

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
  
  /// increment of the plastic deformation gradient over _dt
  /// it is necessary to update the plastic work
  MaterialProperty<RankTwoTensor> & _fp_increment;
  
  /// Plastic work calculated according to
  /// equation (17) in:
  /// Elastic plastic deformation at finite strains
  /// E. H. Lee 1968,
  /// Stanford University technical report AD678483
  MaterialProperty<Real> & _plastic_work;
  
  /// and the value at previous time step
  const MaterialProperty<Real> & _plastic_work_old;
  
  /// prefactor applied to the plastic work to determine the fraction
  /// of plastic energy that contributes to damage.
  const Real _plastic_damage_prefactor;
  
  /// UserObject to read the initial plastic deformation gradient from file
  /// The file will have one row for each element
  /// each row will contain the components
  /// Fp_{11} Fp_{12} Fp_{13} Fp_{21} Fp_{22} Fp_{23} Fp_{31} Fp_{32} Fp_{33}
  /// of the initial plastic deformation gradient
  const ElementPropertyReadFile * const _read_initial_Fp;
  
  /// Coupled temperature
  const VariableValue & _temperature;
  
  /// Reference temperature with no thermal expansion
  const Real _reference_temperature;
  
  /// Linear thermal expansion coefficient
  /// and its first derivative with respect to temperature  
  const Real _thermal_expansion;
  const Real _dCTE_dT;
  
  /// Output old values of pk2 and Fp if NR algorithm fails
  bool _suppress_constitutive_failure;
  
  /// Volumetric thermal expansion
  Real _volumetric_thermal_expansion;
  
  /// Thermal eigenstrain
  RankTwoTensor _thermal_eigenstrain;
  
};
