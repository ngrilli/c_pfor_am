// Nicolò Grilli
// University of Bristol
// 22 Marzo 2022

#include "IterationAdaptiveDTMax.h"

#include <limits>
#include <set>

registerMooseObject("MooseApp", IterationAdaptiveDTMax);

InputParameters
IterationAdaptiveDTMax::validParams()
{
  InputParameters params = IterationAdaptiveDT::validParams();
  params.addClassDescription("Adjust the timestep based on the number of iterations. "
                             "A maximum time step can be set. ");
  params.addParam<Real>("upper_limit_dt",
                        1.0e6,
                        "Upper limit of the time step.");
  return params;
}

IterationAdaptiveDTMax::IterationAdaptiveDTMax(const InputParameters & parameters)
  : IterationAdaptiveDT(parameters),
    _upper_limit_dt(getParam<Real>("upper_limit_dt"))
{
  if (isParamValid("optimal_iterations"))
  {
    _adaptive_timestepping = true;
    _optimal_iterations = getParam<int>("optimal_iterations");

    if (isParamValid("iteration_window"))
      _iteration_window = getParam<int>("iteration_window");
    else
      _iteration_window = ceil(_optimal_iterations / 5.0);
  }
  else
  {
    if (isParamValid("iteration_window"))
      mooseError("'optimal_iterations' must be used for 'iteration_window' to be used");
    if (isParamValid("linear_iteration_ratio"))
      mooseError("'optimal_iterations' must be used for 'linear_iteration_ratio' to be used");
  }

  if (isParamValid("timestep_limiting_function"))
    _max_function_change =
        isParamValid("max_function_change") ? getParam<Real>("max_function_change") : -1;
  else
  {
    if (isParamValid("max_function_change"))
      mooseError("'timestep_limiting_function' must be used for 'max_function_change' to be used");
    if (_force_step_every_function_point)
      mooseError("'timestep_limiting_function' must be used for 'force_step_every_function_point' "
                 "to be used");
  }
}

Real
IterationAdaptiveDTMax::computeDT()
{
  Real dt = _dt_old;

  if (_cutback_occurred)
  {
    _cutback_occurred = false;
    if (_adaptive_timestepping)
    {
      // Don't allow it to grow this step, but shrink if needed
      bool allowToGrow = false;
      computeAdaptiveDT(dt, allowToGrow);
    }
  }
  else if (_tfunc_last_step)
  {
    _tfunc_last_step = false;
    _sync_last_step = false;
    dt = _time_ipol.sample(_time_old);

    if (_verbose)
      _console << "Setting dt to value specified by dt function: " << std::setw(9) << dt
               << std::endl;
  }
  else if (_sync_last_step)
  {
    _sync_last_step = false;
    dt = _dt_old;

    if (_verbose)
      _console << "Setting dt to value used before sync: " << std::setw(9) << dt << std::endl;
  }
  else if (_adaptive_timestepping)
    computeAdaptiveDT(dt);
  else if (_use_time_ipol)
    dt = computeInterpolationDT();
  else
  {
    dt *= _growth_factor;
    if (dt > _dt_old * _growth_factor)
      dt = _dt_old * _growth_factor;
  }
  
  // Limit the maximum value of the time step
  if (dt > _upper_limit_dt) {
    dt = _upper_limit_dt;
  }

  return dt;
}
