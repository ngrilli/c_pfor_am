// Nicolò Grilli
// National University of Singapore
// 7 Marzo 2021

#include "SelfAnnihilation.h"

registerMooseObject("MooseApp", SelfAnnihilation);

template <>
InputParameters
validParams<SelfAnnihilation>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Dislocation annihilation rate. "
                             "Dislocation density annihilates itself.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip velocity, "
							   "for instance from 0 to 11 for FCC.");
  params.addParam<Real>("dc",1.0,"Characteristic length of dislocation annihilation");							   
  return params;
}

SelfAnnihilation::SelfAnnihilation(const InputParameters & parameters)
  : Kernel(parameters),
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dc(getParam<Real>("dc")),
	_dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")) // Velocity value (signed)
{
}

Real
SelfAnnihilation::computeQpResidual()
{
  Real val;
  
  val = std::abs(_dislo_velocity[_qp][_slip_sys_index]);
  val *= _u[_qp];
  val *= _u[_qp];
  val *= _dc;
  
  return _test[_i][_qp] * val;
}

Real
SelfAnnihilation::computeQpJacobian()
{
  Real val;
  
  val = std::abs(_dislo_velocity[_qp][_slip_sys_index]);
  val *= _u[_qp];
  val *= _dc;
  val *= 2.0;
  
  return _test[_i][_qp] * _phi[_j][_qp] * val;
}

// Self annihilation, therefore no off-diagonal term
Real
SelfAnnihilation::computeQpOffDiagJacobian(unsigned int jvar)
{	  
  return 0.0; 
}
