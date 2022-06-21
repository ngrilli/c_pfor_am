// Nicolò Grilli
// University of Bristol
// 4 Giugno 2022

#pragma once

#include "ComputeCrystalPlasticityEigenstrainBase.h"
#include "DerivativeMaterialInterface.h"
#include "ElementPropertyReadFile.h"

/**
 * Compute eigenstrain representing residual deformation
 * that induces residual stress in a sample
 * Eigenstrain grows to the maximum value based on a coupled variable
 * called residual_def_level that grows from 0 to 1
 * The maximum residual deformation can be homogeneous in space
 * or can be read from file to apply a non-homogeneous distribution
 * The residual deformation is represented by 9 numbers
 * because it is a RankTwoTensor
 */
class ComputeCrystalPlasticityResidualEigenstrain
  : public DerivativeMaterialInterface<ComputeCrystalPlasticityEigenstrainBase>
{
public:
  static InputParameters validParams();

  ComputeCrystalPlasticityResidualEigenstrain(const InputParameters & parameters);

  /// We need to set initial values for lattice thermal expansion coefficients
  virtual void initQpStatefulProperties() override;

protected:
  ///Compute the deformation gradient due to thermal expansion
  virtual void computeQpDeformationGradient() override;

  /// Residual deformation level from 0 to 1
  /// Variable value
  const VariableValue & _residual_def_level;

  ///Stores the derivative of the deforamtion gradient w.r.t _residual_def_level
  MaterialProperty<RankTwoTensor> & _ddeformation_gradient_dlevel;

  ///Residual deformation tensor minus identity
  ///in the global reference frame built from the 9 components
  ///provided in the input file, note the matrix is filled in the following order
  ///00, 10, 20, 01, 11, 21, 02, 12, 22
  ///the input file should contain the wanted eigenstrain deformation gradient
  ///minus identity, because identity is added in computeQpDeformationGradient()
  const RankTwoTensor _residual_def_components;

  /**
   * UserObject to read the components of the initial residual deformation from file
   */
  const ElementPropertyReadFile * const _read_initial_residual_def;

  ///Residual deformation tensor in the global reference frame
  ///minus identity matrix
  MaterialProperty<RankTwoTensor> & _residual_def;

};
