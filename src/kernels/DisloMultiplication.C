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
  params.addClassDescription("Dislocation multiplication rate. "
                             "Edge dislocations lead to screw multiplication and vice versa");
  params.addCoupledVar("rho_mult_1", 0.0, "First type of dislocations leading to multiplication of this type.");
  params.addCoupledVar("rho_mult_2", 0.0, "Second type of dislocations leading to multiplication of this type.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip velocity, "
							   "for instance from 0 to 11 for FCC.");
  params.addParam<Real>("Lc",1.0,"Characteristic length of dislocation loops");							   
  return params;
}

DisloMultiplication::DisloMultiplication(const InputParameters & parameters)
  : Kernel(parameters),
    _rho_mult_1(coupledValue("rho_mult_1")),
    _rho_mult_2(coupledValue("rho_mult_2")),
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_Lc(getParam<Real>("Lc")),
    _rho_mult_1_coupled(isCoupled("rho_mult_1")),
	_rho_mult_2_coupled(isCoupled("rho_mult_2")),
    _rho_mult_1_var(_rho_mult_1_coupled ? coupled("rho_mult_1") : 0),
	_rho_mult_2_var(_rho_mult_2_coupled ? coupled("rho_mult_2") : 0),
	_dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")) // Velocity value (signed)
{
}

Real
DisloMultiplication::computeQpResidual()
{
  Real val;
  
  val = std::abs(_dislo_velocity[_qp][_slip_sys_index]);
  val *= (_rho_mult_1[_qp] + _rho_mult_2[_qp]);
  val /= _Lc;
  
  return - _test[_i][_qp] * val;
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
	 
    return - _test[_i][_qp] * _phi[_j][_qp] * std::abs(_dislo_velocity[_qp][_slip_sys_index]) / _Lc;
	
  } else if (_rho_mult_2_coupled && jvar == _rho_mult_2_var) {
	  
    return - _test[_i][_qp] * _phi[_j][_qp] * std::abs(_dislo_velocity[_qp][_slip_sys_index]) / _Lc;
	
  }
  else {
	  
	return 0.0;
	
  } 
}
