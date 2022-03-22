// Nicolò Grilli
// University of Bristol
// 22 Marzo 2022

#pragma once

#include "IterationAdaptiveDT.h"

/**
 * Adjust the timestep based on the number of iterations.
 * The user can specitfy an optimal_iterations number of non-linear iterations and an
 * iteration_window.
 * The time stepper attempts to increase the time step if the non-linear iteration count is below
 * optimal_iterations - iteration_window and it attempts to reduce the time step if the non-linear
 * iteration count is above optimal_iterations + iteration_window. Similar rules apply to the number
 * of linear iterations with the multiplier linear_iteration_ratio applied to optimal_iterations and
 * iteration_window.
 * This time stepper allows the user to specify a function that limits the maximal time step change.
 * This time stepper allows the user to specify a limiting time step length through a postprocessor.
 * A maximum time step can be set.
 */
class IterationAdaptiveDTMax : public IterationAdaptiveDT //public PostprocessorInterface
{
public:
  static InputParameters validParams();

  IterationAdaptiveDTMax(const InputParameters & parameters);

protected:
  virtual Real computeDT() override;

  /// Maximum value of the time step
  const Real & _upper_limit_dt;
  
};

