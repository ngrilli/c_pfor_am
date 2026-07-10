// Daksheena Jeyasingam
// Nicolò Grilli
// Università di Bristol
// 9 Luglio 2026

#include "ComputeMultipleCrystalPlasticityStressNN.h"

registerMooseObject("c_pfor_amApp", ComputeMultipleCrystalPlasticityStressNN);

InputParameters
ComputeMultipleCrystalPlasticityStressNN::validParams()
{
  InputParameters params = ComputeMultipleCrystalPlasticityStress::validParams();

  params.addClassDescription(
      "This class extends ComputeMultipleCrystalPlasticityStress by using a neural network "
      "to provide an initial guess for the Newton iteration. The neural network is trained to predict "
      "the stress increment based on the current state of the material and the deformation gradient. ");
  return params;
}

ComputeMultipleCrystalPlasticityStressNN::ComputeMultipleCrystalPlasticityStressNN(
    const InputParameters & parameters)
  : ComputeMultipleCrystalPlasticityStress(parameters),
    _pk2_initial_output(declareProperty<RankTwoTensor>("pk2_initial")),
    _fp_initial_output(declareProperty<RankTwoTensor>("fp_initial_output")),
    _iteration_output(declareProperty<Real>("iteration_output"))
//    _torch_script_userobject(getUserObject<TorchScriptUserObject>("torch_script_userobject"))
{
  _convergence_failed = false;
}

void
ComputeMultipleCrystalPlasticityStressNN::solveStress()
{
  unsigned int iteration = 0;
  RankTwoTensor dpk2;
  Real rnorm, rnorm0, rnorm_prev;

  // Calculate stress residual
  calculateResidualAndJacobian();
  if (_convergence_failed)
  {
    if (_print_convergence_message)
      mooseWarning("ComputeMultipleCrystalPlasticityStressNN: the slip increment exceeds tolerance at element ",
                   _current_elem->id(),
                   " and Gauss point ",
                   _qp);

    return;
  }

  rnorm = _residual_tensor.L2norm();
  rnorm0 = rnorm;

  // Check for stress residual tolerance; different from user object version which
  // compares the absolute tolerance of only the original rnorm value
  while (rnorm > _rtol * rnorm0 && rnorm > _abs_tol && iteration < _maxiter)
  {
    // Calculate stress increment
    dpk2 = -_jacobian.invSymm() * _residual_tensor;
    _pk2[_qp] = _pk2[_qp] + dpk2;

    _pk2_initial_output[_qp] = _pk2[_qp];

    calculateResidualAndJacobian();

    if (_convergence_failed)
    {
      if (_print_convergence_message)
        mooseWarning("ComputeMultipleCrystalPlasticityStressNN: the slip increment exceeds tolerance "
                     "at element ",
                     _current_elem->id(),
                     " and Gauss point ",
                     _qp);

      return;
    }

    rnorm_prev = rnorm;
    rnorm = _residual_tensor.L2norm();

    if (_use_line_search && rnorm > rnorm_prev && !lineSearchUpdate(rnorm_prev, dpk2))
    {
      if (_print_convergence_message)
        mooseWarning("ComputeMultipleCrystalPlasticityStressNN: failed with line search");

      _convergence_failed = true;
      return;
    }

    if (_use_line_search)
      rnorm = _residual_tensor.L2norm();

    iteration++;
  }

  _iteration_output[_qp] = iteration;

  if (iteration >= _maxiter)
  {
    if (_print_convergence_message)
      mooseWarning("ComputeMultipleCrystalPlasticityStress: Stress Integration error rmax = ",
                   rnorm,
                   " and the tolerance is ",
                   _rtol * rnorm0,
                   " when the rnorm0 value is ",
                   rnorm0,
                   " for element ",
                   _current_elem->id(),
                   " and qp ",
                   _qp);

    _convergence_failed = true;
  }
}
