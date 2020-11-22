// Nicolò Grilli
// National University of Singapore
// 22 Novembre 2020

#include "DisloMultiplication.h"

registerMooseObject("MooseApp", DisloMultiplication);

template <>
InputParameters
validParams<DisloMultiplication>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("rho_mult_1", 0.0, "First type of dislocations leading to multiplication of this type.");
  params.addCoupledVar("rho_mult_2", 0.0, "Second type of dislocations leading to multiplication of this type.");
  return params;
}

DisloMultiplication::DisloMultiplication(const InputParameters & parameters)
  : Kernel(parameters),
    _rho_mult_1(coupledValue("rho_mult_1")),
    _rho_mult_2(coupledValue("rho_mult_2")),
    _rho_mult_1_coupled(isCoupled("rho_mult_1")),
	_rho_mult_2_coupled(isCoupled("rho_mult_2")),
    _rho_mult_1_var(_rho_mult_1_coupled ? coupled("rho_mult_1") : 0),
	_rho_mult_2_var(_rho_mult_2_coupled ? coupled("rho_mult_2") : 0)
{
}

Real
DisloMultiplication::computeQpResidual()
{
  return - _test[_i][_qp] * (_rho_mult_1[_qp] + _rho_mult_2[_qp]);
}

Real
DisloMultiplication::computeQpJacobian()
{
  return 0.0;
}

Real
DisloMultiplication::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_rho_mult_1_coupled && jvar == _rho_mult_1_var)
  {
	 
    return - _test[_i][_qp] * _phi[_j][_qp];
	
  } else if (_rho_mult_2_coupled && jvar == _rho_mult_2_var) {
	  
    return - _test[_i][_qp] * _phi[_j][_qp];
	
  }
  else {
	  
	return 0.0;
	
  } 
}
