// Nicolò Grilli
// University of Bristol
// 4 Giugno 2022

// Compute eigenstrain representing residual deformation
// that induces residual stress in a sample
// Eigenstrain grows to the maximum value based on a coupled variable
// called residual_def_level that grows from 0 to 1 
// The maximum residual deformation can be homogeneous in space
// or can be read from file to apply a non-homogeneous distribution
// The residual deformation is represented by 9 numbers
// because it is a RankTwoTensor

#include "ComputeCrystalPlasticityResidualEigenstrain.h"

registerMooseObject("TensorMechanicsApp", ComputeCrystalPlasticityResidualEigenstrain);

InputParameters
ComputeCrystalPlasticityResidualEigenstrain::validParams()
{
  InputParameters params = ComputeCrystalPlasticityEigenstrainBase::validParams();

  params.addCoupledVar("residual_def_level", "Residual deformation level from 0 to 1");

  // Let's check the range of the parameter here
  params.addRequiredRangeCheckedParam<std::vector<Real>>(
      "residual_def_components",
      "residual_def_components_size=9",
      "Vector of values defining the maximum magnitude of the components" 
	  " of the residual deformation tensor that will be applied when residual_def_level = 1."
	  " Note that these components refers to the beta tensor, meaning the actual eigenstrain"
	  " deformation gradient will be: I+beta."
	  " Note the matrix is filled in the following order:"
	  " 00, 10, 20, 01, 11, 21, 02, 12, 22");

  return params;
}

ComputeCrystalPlasticityResidualEigenstrain::ComputeCrystalPlasticityResidualEigenstrain(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ComputeCrystalPlasticityEigenstrainBase>(parameters),
    _residual_def_level(coupledValue("residual_def_level")),
    _ddeformation_gradient_dlevel(declarePropertyDerivative<RankTwoTensor>(
        _deformation_gradient_name, getVar("residual_def_level", 0)->name())),
    _residual_def_components(getParam<std::vector<Real>>("residual_def_components")),
    _residual_def(declareProperty<RankTwoTensor>(
        _eigenstrain_name +
        "_residual_def")) // avoid duplicated material name by including the eigenstrain name 
{
}

void
ComputeCrystalPlasticityResidualEigenstrain::initQpStatefulProperties()
{
  ComputeCrystalPlasticityEigenstrainBase::initQpStatefulProperties();
  // Assign the constant read residual deformation tensor in the global
  // reference frame to the corresponding material property
  // here identity is added to residual_def_components
  // therefore, the input file should contain the wanted eigenstrain deformation gradient
  // minus identity
  _residual_def[_qp] = _residual_def_components;
}

void
ComputeCrystalPlasticityResidualEigenstrain::computeQpDeformationGradient()
{
  // compute the residual deformation gradient based on _residual_def_level
  // this is homogeneous in the geometry
  // when _residual_def_level[_qp] = 1.0
  // the full eigenstrain is applied, therefore no incremental formulation is needed
  _deformation_gradient[_qp] = RankTwoTensor::Identity() 
                             + _residual_def_level[_qp] * _residual_def[_qp];

  // compute the derivative of eigenstrain deformation gradient w.r.t _residual_def_level
  _ddeformation_gradient_dlevel[_qp] = _residual_def[_qp];
}
