// Nicolò Grilli
// University of Bristol
// 5 Dicembre 2022

#pragma once

#include "Diffusion.h"
#include "Material.h"

// Forward Declarations
class TDepHeatConduction;

/**
 * Computes residual/Jacobian contribution for $(k (T) \\nabla T, \\nabla \\psi)$ term.
 * k(T) depends linearly on temperature. 
 */
class TDepHeatConduction : public Diffusion
{
public:
  static InputParameters validParams();

  TDepHeatConduction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const MaterialProperty<Real> & _diffusion_coefficient;

  /// Constant derivative of the thermal conductivity with respect to temperature.
  const Real _dthermal_conductivity_dT;
  
  /// Reference temperature (K) at which thermal conductivity is _diffusion_coefficient[_qp].
  const Real _reference_temperature;
};
