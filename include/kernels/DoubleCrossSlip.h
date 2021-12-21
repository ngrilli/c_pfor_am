// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// University of Bristol
// 13 Novembre 2021

// Cross slip rate based on equation 2.5 in Christophe Depres 2004 Thesis
// MODÉLISATION PHYSIQUE DES STADES PRÉCURSEURS
// DE L’ENDOMMAGEMENT EN FATIGUE
// DANS L’ACIER INOXYDABLE AUSTÉNITIQUE 316L
// http://www.numodis.fr/download/2004_These_Christophe_DEPRES.pdf

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
  
  // Ratio between Schmid factor of the
  // cross slip system and of the primary system
  const Real _cssf;
  
  // Stage III resolved shear stress
  const Real _tauIII;
  
  // Thermal coefficient of stage III resolved shear stress
  const Real _dtauIII_dT;

  // reference temperature for stage III resolved shear stress
  const Real _reference_temperature;
  
  // Boltzmann constant
  const Real _kB;
  
  // Activation volume
  const Real _Vact;
  
  // Upper limit of the cross slip rate of dislocation density
  // It corresponds to the situation in which most screw dislocations
  // on the primary slip plane would move to a cross slip plane
  const Real _drho_cs_tol;

  // radius of the cross slip dislocation segment
  // it is about 44 burgers vectors according to
  // Nicolò Grilli, Koenraad G.F. Janssens, Helena Van Swygenhoven
  // Crystal plasticity finite element modelling of low cycle fatigue in fcc metals
  // Journal of the Mechanics and Physics of Solids 84 (2015) 424–435
  const Real _R_cs;

  // resolved shear stress value (signed) on all slip systems
  // it is used to calculate the resolved shear stress
  // on the cross slip system
  const MaterialProperty<std::vector<Real>> & _tau_out;

};

