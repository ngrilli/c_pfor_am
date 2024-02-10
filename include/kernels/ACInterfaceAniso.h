// Nicolò Grilli
// Università di Bristol
// 10 Febbraio 2024

#pragma once

#include "ACInterface.h"

/**
 * Compute the anisotropic Allen-Cahn interface term with the weak form residual
 * \f$ \left( \kappa_i \nabla\eta_i, \nabla (L_i \psi) \right) \f$
 * where \kappa_i is calculated based on crystal orientation, following equations (12) and (16) in:
 * Min Yang, Lu Wang and Wentao Yan
 * Phase-ﬁeld modeling of grain evolutions in additive
 * manufacturing from nucleation, growth, to coarsening
 * npj Computational Materials volume 7, 56 (2021)
 */
class ACInterfaceAniso : public ACInterface
{
public:
  static InputParameters validParams();

  ACInterfaceAniso(const InputParameters & parameters);

protected:
  virtual Real computeAnisotropy();

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  
  /// Grain boundary energy anisotropy coefficient
  const Real _e_anisotropy;

};
