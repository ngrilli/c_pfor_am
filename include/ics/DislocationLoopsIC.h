// Nicolò Grilli
// 1 Febbraio 2021
// National University of Singapore

#pragma once

#include "InitialCondition.h"

// Forward Declarations
class DislocationLoopsIC;
class InputParameters;

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<DislocationLoopsIC>();

/**
 * Initialize dislocation loops.
 */
class DislocationLoopsIC : public InitialCondition
{
public:
  static InputParameters validParams();

  DislocationLoopsIC(const InputParameters & parameters);

protected:

  virtual Real value(const Point & p) override;
  
  std::vector<Real> _centrex;
  std::vector<Real> _centrey;
  std::vector<Real> _radii;
  std::vector<Real> _width;
  std::vector<Real> _rho_max;
  
  // Type of variable for the dislocation loop
  MooseEnum _variable_type;

};
