// Nicolò Grilli
// Università di Bristol
// 18 Aprile 2026

#pragma once

#include "Material.h"

/**
 * A material property that uses a piecewise linear function of a coupled variable
 */
class VariablePiecewiseLinearMaterial : public Material
{
public:
  static InputParameters validParams();

  VariablePiecewiseLinearMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Coupled variable
  const VariableValue & _var;

  /// x-y tabular data
  std::vector<Real> _x;
  std::vector<Real> _y;

  /// The material property computed by this material, which is a piecewise linear function of the coupled variable
  MaterialProperty<Real> & _f;
};
