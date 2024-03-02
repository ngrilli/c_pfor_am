// Nicolò Grilli
// Parsa Esmati
// Università di Bristol
// Fernando Valiente Dies
// ANSTO
// 10 Dicembre 2023

#pragma once

#include "InputParameters.h"
#include "Action.h"

/**
 * Action that sets up TimeDerivative, Reaction, BodyForce, CoupledTanh, Diffusion, CoupledPhaseGrain kernels,
 * for the zeta variable, which is 0 in the liquid phase and 1 in the solid phase.
 */
class LiquidSolidKernelAction : public Action
{
public:
  static InputParameters validParams();

  LiquidSolidKernelAction(const InputParameters & params);

  virtual void act();

protected:
  
  /// number of grains
  const unsigned int _op_num;

  /// base name for the order parameter variables
  const std::string _var_name_base;
  
  /// Model parameters
  const Real _sigma_p;
  const Real _delta_f_p;
  const Real _l_p;
  const Real _sigma_g0;
  const Real _delta_f_g;
  const Real _l_g;
  const Real _L_p;
  const Real _theta;
  const Real _T_l;
  
  /// Interaction coefficient between zeta and eta variables
  const Real _gamma_p;
  
  Real _m_p;
  Real _m_g;
  Real _k_p;
  Real _k_g;
  
};
