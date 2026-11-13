//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMVectorComponentAux.h"

registerMooseObject("c_pfor_amApp", MFEMVectorComponentAux);

InputParameters
MFEMVectorComponentAux::validParams()
{
  InputParameters params = MFEMAuxKernel::validParams();
  params.addClassDescription("Projects a single component of a vector coefficient onto a scalar MFEM auxvariable.");
  params.addRequiredParam<MFEMVectorCoefficientName>("vector_coefficient","Name of the vector coefficient to extract a component from.");
  params.addRequiredParam<int>("component","Zero-based component index to extract (0, 1, 2 for x, y, z in 3D).");
  return params;
}

MFEMVectorComponentAux::MFEMVectorComponentAux(const InputParameters & parameters)
  : MFEMAuxKernel(parameters),
    _vec_coef(getVectorCoefficient("vector_coefficient")),
    _component(getParam<int>("component")),
    _component_coef(_vec_coef, _component)
{
}

void
MFEMVectorComponentAux::execute()
{
  _result_var.ProjectCoefficient(_component_coef);
}

#endif