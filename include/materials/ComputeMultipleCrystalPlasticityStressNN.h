// Daksheena Jeyasingam
// Nicolò Grilli
// Università di Bristol
// 9 Luglio 2026

#pragma once

#include "ComputeMultipleCrystalPlasticityStress.h"
#include "TorchScriptUserObject.h"

/**
 * This class extends ComputeMultipleCrystalPlasticityStress by using a neural network
 * to provide an initial guess for the Newton iteration. The neural network is trained to predict
 * the stress increment based on the current state of the material and the deformation gradient. 
 */
class ComputeMultipleCrystalPlasticityStressNN : public ComputeMultipleCrystalPlasticityStress
{
public:
  static InputParameters validParams();

  ComputeMultipleCrystalPlasticityStressNN(const InputParameters & parameters);

protected:
  void solveStress();

  // Material properties exposed for AuxKernels / Exodus / CSV output.
  MaterialProperty<RankTwoTensor> & _pk2_initial_output;
  MaterialProperty<RankTwoTensor> & _fp_initial_output;
  MaterialProperty<Real> & _iteration_output;

  /// The user object that holds the torch module
//  const TorchScriptUserObject & _torch_script_userobject;
};
