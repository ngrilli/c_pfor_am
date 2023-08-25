// Nicolò Grilli
// Università di Bristol
// 25 Agosto 2023

#pragma once

#include "ACGrGrBase.h"

// Forward Declarations

/**
 * Term representing the interaction between phases and grain orientations
 * last term in equation (10) in:
 * Min Yang, Lu Wang, Wentao Yan
 * Phase-field modeling of grain evolutions in additive manufacturing from nucleation, growth, to coarsening
 * https://www.nature.com/articles/s41524-021-00524-6
 */
class GrainSolidification : public ACGrGrBase
{
public:
  static InputParameters validParams();

  GrainSolidification(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _gamma;
  
  /// Phase field representing liquid (0) or solid (1)
  const VariableValue & _zeta;
  const bool _zeta_coupled;
  unsigned int _zeta_var;
};
