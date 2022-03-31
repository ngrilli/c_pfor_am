// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// University of Bristol
// 31 Marzo 2022

#pragma once

// MOOSE includes
#include "AuxKernel.h"

// Calculate directional derivative along edge and screw dislocation propagation direction.
// This AuxKernel applies to a vector auxiliary variable

class ArrayDirectionalDerivative : public ArrayAuxKernel
{
public:
  static InputParameters validParams();

  /**
   * Class constructor
   * @param parameters Input parameters for the object
   */
  ArrayDirectionalDerivative(const InputParameters & parameters);

protected:
  virtual RealEigenVector computeValue() override;

private:

  // The vector variable from which to compute the directional derivative
  const ArrayVariableGradient & _grad_variable;
  
  // Edge slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _edge_slip_direction;

  // Screw slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
  // Character of dislocations (edge or screw)
  const enum class DisloCharacter { edge, screw } _dislo_character;
  
  // advection velocity
  std::vector<Real> _velocity;

};
