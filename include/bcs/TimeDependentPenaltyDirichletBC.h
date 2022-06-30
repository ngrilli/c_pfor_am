// Nicolò Grilli
// University of Bristol
// 30 Giugno 2022

#pragma once

#include "IntegratedBC.h"

class Function;

/**
 * A different approach to applying Dirichlet BCs
 *
 * uses \f$ \int(p u \cdot \phi)=\int(p f \cdot \phi)\f$ on \f$d\omega\f$
 *
 * The penalty coefficient can also be a time and space-dependent function.
 * It is useful to remove Dirichlet BC at specific times,
 * for instance after stress relaxation.
 */

class TimeDependentPenaltyDirichletBC : public IntegratedBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  static InputParameters validParams();

  TimeDependentPenaltyDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const Function & _func;
  const Function & _p_func;
};
