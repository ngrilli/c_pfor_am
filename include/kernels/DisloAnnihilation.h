// Nicolò Grilli
// National University of Singapore
// 25 Novembre 2020

#pragma once

#include "Kernel.h"

class DisloAnnihilation;

template <>
InputParameters validParams<DisloAnnihilation>();

class DisloAnnihilation : public Kernel
{
public:
  static InputParameters validParams();

  DisloAnnihilation(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  // Opposite signed dislocation type leading to annihilation
  const VariableValue & _rho_annih;
  
  // Slip system index to determine slip system velocity
  const unsigned int _slip_sys_index;
  
  // Characteristic length of dislocation loops
  const Real _dc;
  
  const bool _rho_annih_coupled;
  
  unsigned int _rho_annih_var;
  
  // Dislocation velocity value (signed) on all slip systems
  const MaterialProperty<std::vector<Real>> & _dislo_velocity;

};

