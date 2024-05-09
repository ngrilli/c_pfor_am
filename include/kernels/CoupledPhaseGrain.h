// Fernando Valiente Dies
// ANSTO
// Nicolò Grilli
// Parsa Esmati
// Università di Bristol
// 17 Dicembre 2023

#pragma once

#include "Kernel.h"

// Interaction between phases and grain orientations for zeta term.

class CoupledPhaseGrain : public Kernel
{
public:
  static InputParameters validParams();

  CoupledPhaseGrain(const InputParameters & parameters);

protected:
  virtual std::vector<Real> assignOps();

  virtual Real computeQpResidual() override;
  
  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// The coupled order parameters
  const unsigned int _op_num;

  const std::vector<const VariableValue *> _vals;
  const std::vector<unsigned int> _vals_var;
  
  /// Kinetic coefficient for phase transformation
  const Real _L_p;
  
  /// Interaction coefficient between zeta and eta variables
  const Real _gamma_p;
  
  /// The m_g parameter of the free energy calculated in GBEvolution
  const MaterialProperty<Real> & _m_g;

};
