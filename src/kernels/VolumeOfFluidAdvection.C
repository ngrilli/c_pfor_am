// Nicolò Grilli
// University of Bristol
// 29 Settembre 2022

#include "VolumeOfFluidAdvection.h"

registerMooseObject("LevelSetApp", VolumeOfFluidAdvection);

InputParameters
VolumeOfFluidAdvection::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addClassDescription("Implements the level set advection equation: $\\vec{v}\\cdot\\nabla "
                             "u = 0$, where the weak form is $(\\psi_i, \\vec{v}\\cdot\\nabla u) = "
                             "0$. The difference with LevelSetAdvection is the usage of three coupled "
                             "variables for the components of the velocity.");
  params.addRequiredCoupledVar("u", "The x velocity variable.");
  params.addRequiredCoupledVar("v", "The y velocity variable.");
  params.addCoupledVar("w", "The z velocity variable.");
  return params;
}

VolumeOfFluidAdvection::VolumeOfFluidAdvection(const InputParameters & parameters)
  : ADKernelValue(parameters), 
    _u(coupledValue("u")),
    _v(coupledValue("v")),
    _has_w(isCoupled("w")),
    _w(_has_w ? coupledValue("w") : _zero)
{
}

ADReal
VolumeOfFluidAdvection::precomputeQpResidual()
{
  RealVectorValue velocity;
  
  velocity(0) = _u[_qp];
  velocity(1) = _v[_qp];
  
  if (_has_w) {
	  
    velocity(2) = _w[_qp];
	  
  } else {
	  
    velocity(2) = 0.0;
	  
  }
	
  return velocity * _grad_u[_qp];
}
