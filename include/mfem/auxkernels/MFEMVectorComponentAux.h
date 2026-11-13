//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "MFEMAuxKernel.h"

/**
 * Projects a single component of a vector coefficient onto a scalar MFEM auxvariable.
 *
 * The component index is zero-based (0, 1, 2 for x, y, z components in 3D).
 */
class MFEMVectorComponentAux : public MFEMAuxKernel
{
public:
  static InputParameters validParams();

  MFEMVectorComponentAux(const InputParameters & parameters);

  virtual ~MFEMVectorComponentAux() = default;

  virtual void execute() override;

protected:
  /// Reference to source vector coefficient.
  mfem::VectorCoefficient & _vec_coef;

  /// Component index to extract (0-based).
  const int _component;

  /// Coefficient that extracts the requested component from the vector coefficient.
  mfem::VectorComponentCoefficient _component_coef;
};

#endif