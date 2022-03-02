// Nicolò Grilli
// University of Bristol
// 2 Luglio 2021

// DG upwinding for the advection of a coupled variable
// It is assumed that upwind scheme is always used
// Therefore, no need for a upwinding_type flag as in ConservativeAdvectionCoupled
// This kernel implements the term (on LHS): 
// d(rho_coupled v)/dx if dislo_character = edge
// d(rho_coupled v)/dy if dislo_character = screw
// with constant, user-defined, uniform velocity field

#pragma once

#include "DGKernel.h"

class DGAdvectionCoupledVConst : public DGKernel
{
public:
  static InputParameters validParams();

  DGAdvectionCoupledVConst(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;
  
  // Coupled dislocation density in the flux term
  const VariableValue & _rho_coupled;
  
  const bool _rho_coupled_coupled;
  unsigned int _rho_coupled_var;
  
  // Coupled dislocation density in the neighbouring element
  const VariableValue & _rho_neighbor;

  /// advection velocity
  RealVectorValue _velocity;

};

