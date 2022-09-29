// Nicolò Grilli
// University of Bristol
// 29 Settembre 2022

#include "VoFTimeDerivativeSUPG.h"

registerMooseObject("LevelSetApp", VoFTimeDerivativeSUPG);

InputParameters
VoFTimeDerivativeSUPG::validParams()
{
  InputParameters params = ADTimeKernelGrad::validParams();
  params.addClassDescription(
      "SUPG stablization terms for the time derivative of the level set equation. "
      "The difference with LevelSet TimeDerivativeSUPG is the usage of three coupled "
      "variable for the components of the velocity. ");
  params.addRequiredCoupledVar("u", "The x velocity variable.");
  params.addRequiredCoupledVar("v", "The y velocity variable.");
  params.addCoupledVar("w", "The z velocity variable.");
  return params;
}

VoFTimeDerivativeSUPG::VoFTimeDerivativeSUPG(const InputParameters & parameters)
  : ADTimeKernelGrad(parameters), 
    _u(coupledValue("u")),
    _v(coupledValue("v")),
    _has_w(isCoupled("w")),
    _w(_has_w ? coupledValue("w") : _zero)
{
}

ADRealVectorValue
VoFTimeDerivativeSUPG::precomputeQpResidual()
{
  RealVectorValue velocity;
  
  velocity(0) = _u[_qp];
  velocity(1) = _v[_qp];
  
  if (_has_w) {
	  
    velocity(2) = _w[_qp];
	  
  } else {
	  
    velocity(2) = 0.0;
	  
  }	
	
  ADReal tau =
      _current_elem->hmin() /
      (2 * (velocity + RealVectorValue(libMesh::TOLERANCE * libMesh::TOLERANCE)).norm());
  return tau * velocity * _u_dot[_qp];
}
