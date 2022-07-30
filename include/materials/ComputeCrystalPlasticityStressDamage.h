// Nicolò Grilli
// University of Bristol
// 17 Luglio 2022

#pragma once

#include "ComputeFiniteStrainElasticStress.h"

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
 *
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
class ComputeCrystalPlasticityStressDamage : public ComputeFiniteStrainElasticStress
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
  virtual void updateStress(RankTwoTensor & cauchy_stress, RankFourTensor & jacobian_mult);

  /**
   * initializes the stateful properties such as PK2 stress, resolved shear
   * stress, plastic deformation gradient, slip system resistances, etc.
   * This class is often overwritten by inherting classes.
   */
  virtual void initQpStatefulProperties() override;

  /**
   * Calls the residual and jacobian functions used in the stress update
   * algorithm.
   */
  void calculateResidualAndJacobian();

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
   * Solves the internal variables stress as a function of the slip specified
   * by the constitutive model defined in the inheriting class
   */
  void solveStateVariables();

  /**
   * solves for stress, updates plastic deformation gradient.
   */
  void solveStress();

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
									   
  // compute history variable and assign to _E
  // which is used by the fracture model for damage growth
  // Damage grows only because of the positive part of the elastic energy F_pos									 
  virtual void computeHistoryVariable(Real & F_pos, Real & F_neg);

  ///@{Calculates the tangent moduli for use as a preconditioner, using the elastic or elastic-plastic option as specified by the user
  void calcTangentModuli(RankFourTensor & jacobian_mult);
  void elasticTangentModuli(RankFourTensor & jacobian_mult);
  void elastoPlasticTangentModuli(RankFourTensor & jacobian_mult);
  ///@}

  /// performs the line search update
  bool lineSearchUpdate(const Real & rnorm_prev, const RankTwoTensor & dpk2);

  /**
   * Calculates the deformation gradient due to eigenstrain
   */
  void calculateEigenstrainDeformationGrad();

  /// number of plastic models
  const unsigned _num_models;

  /// The user supplied cyrstal plasticity consititutive models
  // CrystalPlasticityDislocationUpdateBase is used in which
  // the method calculateSchmidTensor is virtual
  std::vector<CrystalPlasticityDislocationUpdateBase *> _models;

  /// number of eigenstrains
  const unsigned _num_eigenstrains;

  /// The user supplied cyrstal plasticity eigenstrains
  std::vector<ComputeCrystalPlasticityEigenstrainBase *> _eigenstrains;
  
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
  
  // Plastic work calculated according to
  // equation (17) in:
  // Elastic plastic deformation at finite strains
  // E. H. Lee 1968,
  // Stanford University technical report AD678483
  MaterialProperty<Real> & _plastic_work;
  
  // and the value at previous time step
  const MaterialProperty<Real> & _plastic_work_old;
  
  /// Elastic deformation gradient needed by the user object
  /// PlasticWorkRate to calculate the plastic work
  MaterialProperty<RankTwoTensor> & _elastic_deformation_grad;
  
  // prefactor applied to the plastic work to determine the fraction
  // of plastic energy that contributes to damage.
  const Real _plastic_damage_prefactor;

  /// optional parameter to define several mechanical systems on the same block, e.g. multiple phases
  const std::string _base_name;

  /// Elasticity tensor as defined by a separate class
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  /// Stress residual equation relative tolerance
  Real _rtol;
  /// Stress residual equation absolute tolerance
  Real _abs_tol;

  /// Residual tensor
  RankTwoTensor _residual_tensor;
  /// Jacobian tensor
  RankFourTensor _jacobian;

  /// Maximum number of iterations for stress update
  unsigned int _maxiter;
  /// Maximum number of iterations for internal variable update
  unsigned int _maxiterg;

  /// Type of tangent moduli calculation
  const enum class TangentModuliType { EXACT, NONE } _tan_mod_type;

  /// Maximum number of substep iterations
  unsigned int _max_substep_iter;

  /// time step size during substepping
  Real _substep_dt;

  /// Flag to activate line serach
  bool _use_line_search;

  /// Minimum line search step size
  Real _min_line_search_step_size;

  /// Line search bisection method tolerance
  Real _line_search_tolerance;

  /// Line search bisection method maximum iteration number
  unsigned int _line_search_max_iterations;

  /// strain formulation
  const enum class LineSearchMethod { CutHalf, Bisection } _line_search_method;

  ///@{Plastic deformation gradient RankTwoTensor for the crystal
  MaterialProperty<RankTwoTensor> & _plastic_deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _plastic_deformation_gradient_old;
  ///@}

  ///@{ Generalized eigenstrain deformation gradient RankTwoTensor for the crystal
  MaterialProperty<RankTwoTensor> * _eigenstrain_deformation_gradient;
  const MaterialProperty<RankTwoTensor> * _eigenstrain_deformation_gradient_old;
  ///@}

  ///@{Total deformation gradient RankTwoTensor for the crystal
  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  ///@}

  ///@{Second Piola-Kirchoff stress measure
  MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<RankTwoTensor> & _pk2_old;
  ///@}

  /// Lagrangian total strain measure for the entire crystal
  MaterialProperty<RankTwoTensor> & _total_lagrangian_strain;

  /**
   * Tracks the rotation of the crystal during deformation
   * Note: this rotation tensor is not applied to the crystal lattice
   */
  MaterialProperty<RankTwoTensor> & _updated_rotation;

  /**
   * Crystal rotation in the original, or reference, configuration as defined by
   * Euler angle arguments in the ComputeElasticityTensor classes
   */
  const MaterialProperty<RankTwoTensor> & _crysrot;

  ///@{Helper deformation gradient tensor variables used in iterative solve
  RankTwoTensor _temporary_deformation_gradient;
  RankTwoTensor _elastic_deformation_gradient;
  RankTwoTensor _inverse_plastic_deformation_grad;
  RankTwoTensor _inverse_plastic_deformation_grad_old;
  RankTwoTensor _inverse_eigenstrain_deformation_grad;
  ///@}

  /// Flag to print to console warning messages on stress, constitutive model convergence
  const bool _print_convergence_message;
  
  // UserObject to read the initial plastic deformation gradient from file
  // The file will have one row for each element
  // each row will contain the components
  // Fp_{11} Fp_{12} Fp_{13} Fp_{21} Fp_{22} Fp_{23} Fp_{31} Fp_{32} Fp_{33}
  // of the initial plastic deformation gradient
  const ElementPropertyReadFile * const _read_initial_Fp;

  /// Flag to check whether convergence is achieved or if substepping is needed
  bool _convergence_failed;

  ///@{ Used for substepping; Uniformly divides the increment in deformation gradient
  RankTwoTensor _delta_deformation_gradient;
  RankTwoTensor _temporary_deformation_gradient_old;
  ///@}

  /// Scales the substepping increment to obtain deformation gradient at a substep iteration
  Real _dfgrd_scale_factor;
};
