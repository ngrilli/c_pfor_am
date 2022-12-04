// Nicolò Grilli
// University of Bristol
// 4 Dicembre 2022

#pragma once

// MOOSE includes
#include "TimeDerivative.h"
#include "Material.h"

// Forward Declarations

/**
 * A class for defining the time derivative of the heat equation.
 *
 * By default this Kernel computes:
 *   \f$ \rho * c_p (T) * \frac{\partial T}{\partial t}, \f$
 * where \f$ \rho \f$ and \f$ c_p \f$ are material properties with the names "density" and
 * "specific_heat", respectively.
 * The specific heat depends linearly on temperature
 */
class TDepCpHeatTimeDerivative : public TimeDerivative
{
public:
  /// Contructor for Heat Equation time derivative term.
  static InputParameters validParams();

  TDepCpHeatTimeDerivative(const InputParameters & parameters);

protected:
  /// Compute the residual of the Heat Equation time derivative.
  virtual Real computeQpResidual();

  /// Compute the jacobian of the Heat Equation time derivative.
  virtual Real computeQpJacobian();

  const MaterialProperty<Real> & _specific_heat;
  const MaterialProperty<Real> & _density;
  
  /// Constant derivative of the specific heat with respect to temperature.
  const Real _dspecific_heat_dT;
  
  /// Reference temperature (K) at which specific heat is _specific_heat[_qp]
  const Real _reference_temperature;
};
