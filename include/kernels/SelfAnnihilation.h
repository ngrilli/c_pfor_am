// Nicolò Grilli
// National University of Singapore
// 7 Marzo 2021

#pragma once

#include "Kernel.h"

class SelfAnnihilation;

template <>
InputParameters validParams<SelfAnnihilation>();

class SelfAnnihilation : public Kernel
{
public:
  static InputParameters validParams();

  SelfAnnihilation(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  
  // Slip system index to determine slip system velocity
  const unsigned int _slip_sys_index;
  
  // Characteristic length of dislocation loops
  const Real _dc;

  // Dislocation velocity value (signed) on all slip systems
  const MaterialProperty<std::vector<Real>> & _dislo_velocity;

};

