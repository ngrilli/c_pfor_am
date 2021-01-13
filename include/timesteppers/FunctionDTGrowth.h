// Nicolò Grilli
// 13 Gennaio 2021
// National University of Singapore

#pragma once

#include "TimeStepper.h"
#include "FunctionInterface.h"
#include "LinearInterpolation.h"

class FunctionDTGrowth;
class Function;

template <>
InputParameters validParams<FunctionDTGrowth>();

class FunctionDTGrowth : public TimeStepper, public FunctionInterface
{
public:
  static InputParameters validParams();

  FunctionDTGrowth(const InputParameters & parameters);

  virtual void init() override;

  virtual void postStep() override;
  virtual void rejectStep() override;

protected:
  virtual Real computeInitialDT() override;
  virtual Real computeDT() override;

  void removeOldKnots();

  const std::vector<Real> & _time_t;
  const std::vector<Real> & _time_dt;

  /// true, if we are using `_function`, false if we are using _time_ipol
  bool _use_function;
  /// The time-dependent function specifying the time step size (turn this into a reference then
  /// `time_t` and `time_dt` is removed)
  const Function * _function;

  /// Piecewise linear definition of time stepping
  std::unique_ptr<LinearInterpolation> _time_ipol;

  Real _growth_factor;
  /// True if cut back of the time step occurred
  bool _cutback_occurred;
  Real _min_dt;

  /// Whether or not to interpolate DT between times
  bool _interpolate;

  std::vector<Real> _time_knots;

};
