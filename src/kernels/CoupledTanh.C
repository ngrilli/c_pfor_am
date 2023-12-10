// Nicolò Grilli
// Università di Bristol
// 10 Dicembre 2023

#include "CoupledTanh.h"

registerMooseObject("c_pfor_amApp", CoupledTanh);

InputParameters
CoupledTanh::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Hyperbolic tangent of a coupled variable. ");
  params.addCoupledVar("v", 0.0, "The coupled variable. ");
  params.addParam<Real>("A",0.0,"Prefactor. ");
  params.addParam<Real>("theta",0.0,"Prefactor inside tanh. ");
  params.addParam<Real>("vn",1.0,"Normalization constant for the coupled variable. ");
  return params;
}

CoupledTanh::CoupledTanh(const InputParameters & parameters)
  : Kernel(parameters), 
  _v(coupledValue("v")),
  _v_coupled(isCoupled("v")),
  _v_var(_v_coupled ? coupled("v") : 0),
  _A(getParam<Real>("A")),
  _theta(getParam<Real>("theta")),
  _vn(getParam<Real>("vn"))
{
}

Real
CoupledTanh::computeQpResidual()
{
  return _test[_i][_qp] * _A * std::tanh(_theta * (_v[_qp]/_vn - 1.0));
}

// Variable _u does not appear in residual, therefore 0 on-diagonal Jacobian
Real
CoupledTanh::computeQpJacobian()
{
  return 0.0;
}

Real
CoupledTanh::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real jac;	
	
  if (_v_coupled && jvar == _v_var)
  {
	  
    jac = _test[_i][_qp] * _phi[_j][_qp] *
	       (_A * _theta / _vn) / 
	       std::pow(std::cosh(_theta * (_v[_qp]/_vn - 1.0)),2.0);

    return jac;
	
  }
  else {
	  
	return 0.0;
	
  } 	
}


