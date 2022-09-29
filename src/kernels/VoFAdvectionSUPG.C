// Nicolò Grilli
// University of Bristol
// 29 Settembre 2022

#include "VoFAdvectionSUPG.h"

registerMooseObject("LevelSetApp", VoFAdvectionSUPG);

InputParameters
VoFAdvectionSUPG::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addClassDescription(
      "SUPG stablization term for the advection portion of the level set equation. "
      "The difference with LevelSet TimeDerivativeSUPG is the usage of three coupled "
      "variable for the components of the velocity. ");
  params.addRequiredCoupledVar("u", "The x velocity variable.");
  params.addRequiredCoupledVar("v", "The y velocity variable.");
  params.addCoupledVar("w", "The z velocity variable.");
  return params;
}

VoFAdvectionSUPG::VoFAdvectionSUPG(const InputParameters & parameters)
  : ADKernelGrad(parameters),
    _u(coupledValue("u")),
    _v(coupledValue("v")),
    _has_w(isCoupled("w")),
    _w(_has_w ? coupledValue("w") : _zero)
{
}

ADRealVectorValue
VoFAdvectionSUPG::precomputeQpResidual()
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
  return (tau * velocity) * (velocity * _grad_u[_qp]);
}
