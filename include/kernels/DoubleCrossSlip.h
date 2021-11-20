// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// University of Bristol
// 13 Novembre 2021

#pragma once

#include "Kernel.h"

class DoubleCrossSlip;

template <>
InputParameters validParams<DoubleCrossSlip>();

class DoubleCrossSlip : public Kernel
{
public:
  static InputParameters validParams();

  DoubleCrossSlip(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
private:
  // Slip system index to determine slip system velocity
  const unsigned int _slip_sys_index;
  
  // Edge GND dislocation density  
  const VariableValue & _rho_gnd_edge;

  const bool _rho_gnd_edge_coupled;
  unsigned int _rho_gnd_edge_var;

  // Screw GND dislocation density
  const VariableValue & _rho_gnd_screw;
  
  const bool _rho_gnd_screw_coupled;
  unsigned int _rho_gnd_screw_var;
  
  // Total dislocation density: rho_t 
  const VariableValue & _rho_tot;
  
  const bool _rho_tot_coupled;
  unsigned int _rho_tot_var;

  // Temperature: needed because this process is thermally activated
  const VariableValue & _temp;
  
  // Probability rate for cross slip in El-Azab 2016 paper
  // This is the pre-factor that multiplies total screw
  // density and gives rate of increase of curvature density
  const Real _p_cs;
  
  // Tolerance on small values of remain_rho_tot
  const Real _remain_rho_tol;
  
  // Dislocation mobility
  // used to calculate the resolved shear stress
  const Real _dislo_mobility;
  
  // Ratio between Schmid factor of the
  // cross slip system and of the primary system
  const Real _cssf;
  
  // Stage III resolved shear stress
  const Real _tauIII;
  
  // Boltzmann constant
  const Real _kB;
  
  // Activation volume
  const Real _Vact;

  // Dislocation velocity value (signed) on all slip systems
  const MaterialProperty<std::vector<Real>> & _dislo_velocity;

};

