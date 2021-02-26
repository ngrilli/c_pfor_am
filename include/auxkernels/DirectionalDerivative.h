// Nicolo Grilli
// National University of Singapore
// 26 Febbraio 2021

#pragma once

// MOOSE includes
#include "AuxKernel.h"

// Forward declarations
class DirectionalDerivative;

template <>
InputParameters validParams<DirectionalDerivative>();

/**
 * Extract a component from the gradient of a variable
 */
class DirectionalDerivative : public AuxKernel
{
public:
  static InputParameters validParams();

  /**
   * Class constructor
   * @param parameters Input parameters for the object
   */
  DirectionalDerivative(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:
  /// Reference to the gradient of the coupled variable
  const VariableGradient & _gradient;
  
  // Slip system index to determine slip direction
  const unsigned int _slip_sys_index;
  
  // Edge slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _edge_slip_direction;

  // Screw slip directions of all slip systems
  const MaterialProperty<std::vector<Real>> & _screw_slip_direction;
  
  // Character of dislocations (edge or screw)
  const enum class DisloCharacter { edge, screw } _dislo_character;
  
  /// advection velocity
  std::vector<Real> _velocity;

};

