// Nicolò Grilli
// Università di Bristol
// 10 Febbraio 2024

#pragma once

#include "ACInterface.h"
#include "DelimitedFileReader.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

/**
 * Compute the anisotropic Allen-Cahn interface term with the weak form residual
 * \f$ \left( \kappa_i \nabla\eta_i, \nabla (L_i \psi) \right) \f$
 * where \kappa_i is calculated based on crystal orientation, following equations (12) and (16) in:
 * Min Yang, Lu Wang and Wentao Yan
 * Phase-ﬁeld modeling of grain evolutions in additive
 * manufacturing from nucleation, growth, to coarsening
 * npj Computational Materials volume 7, 56 (2021)
 */
class ACInterfaceAniso : public ACInterface
{
public:
  static InputParameters validParams();

  ACInterfaceAniso(const InputParameters & parameters);

protected:

  virtual void assignEulerAngles();
  
  virtual void computeRotationMatrix();

  virtual Real computeAnisotropy();

  virtual Real computeQpResidual();
  
  /// Grain boundary energy anisotropy coefficient
  const Real _e_anisotropy;
  
  /// Phase field index used to select the Euler angles
  /// and total number of phase fields
  const unsigned int _op;
  const unsigned int _op_num;
  
  /// Euler angles file (in degrees) 
  /// each row must contain three Euler angles 
  /// which correspond to each grain orientation
  std::string _Euler_angles_file_name;
  
  /// Flag to activate continuous model for anisotropy
  /// which does not have a threshold on the gradient of eta_i
  /// for the identification of grain boundaries
  const bool _continuous_anisotropy;
  
  /// Euler angles in degrees for the current grain orientation
  /// and rotation matrices
  RealVectorValue _Euler_angles;
  RotationTensor _R;
  RankTwoTensor _crysrot;

};
