// Nicolò Grilli
// University of Bristol
// 30 Giugno 2022

#include "TimeDependentPenaltyDirichletBC.h"
#include "Function.h"

registerMooseObject("MooseApp", TimeDependentPenaltyDirichletBC);

InputParameters
TimeDependentPenaltyDirichletBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription(
      "Enforces a (possibly) time and space-dependent MOOSE Function Dirichlet boundary condition "
      "in a weak sense by penalizing differences between the current "
      "solution and the Dirichlet data. "
      "The penalty coefficient can also be a time and space-dependent function. "
      "It is useful to remove Dirichlet BC at specific times, "
      "for instance after stress relaxation. ");
  params.addRequiredParam<FunctionName>("forcing_function", "Forcing function");
  params.addRequiredParam<FunctionName>("penalty_function", "Penalty function");

  return params;
}

TimeDependentPenaltyDirichletBC::TimeDependentPenaltyDirichletBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
  _func(getFunction("forcing_function")),
  _p_func(getFunction("penalty_function"))
{
}

Real
TimeDependentPenaltyDirichletBC::computeQpResidual()
{
  return _p_func.value(_t, _q_point[_qp]) * _test[_i][_qp] * (-_func.value(_t, _q_point[_qp]) + _u[_qp]);
}

Real
TimeDependentPenaltyDirichletBC::computeQpJacobian()
{
  return _p_func.value(_t, _q_point[_qp]) * _phi[_j][_qp] * _test[_i][_qp];
}
