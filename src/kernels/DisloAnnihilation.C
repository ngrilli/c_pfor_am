// Nicolò Grilli
// National University of Singapore
// 25 Novembre 2020

#include "DisloAnnihilation.h"

registerMooseObject("MooseApp", DisloAnnihilation);

template <>
InputParameters
validParams<DisloAnnihilation>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Dislocation annihilation rate. "
                             "Positive dislocations annihilate negative ones");
  params.addCoupledVar("rho_annih", 0.0, "Opposite signed dislocation type leading to annihilation.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip velocity, "
							   "for instance from 0 to 11 for FCC.");
  params.addParam<Real>("dc",1.0,"Characteristic length of dislocation annihilation");							   
  return params;
}

DisloAnnihilation::DisloAnnihilation(const InputParameters & parameters)
  : Kernel(parameters),
    _rho_annih(coupledValue("rho_annih")),
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dc(getParam<Real>("dc")),
    _rho_annih_coupled(isCoupled("rho_annih")),
    _rho_annih_var(_rho_annih_coupled ? coupled("rho_annih") : 0),
	_dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")) // Velocity value (signed)
{
}

Real
DisloAnnihilation::computeQpResidual()
{
  Real val;
  
  val = std::abs(_dislo_velocity[_qp][_slip_sys_index]);
  val *= _rho_annih[_qp];
  val *= _u[_qp];
  val *= _dc;
  val *= 4.0;
  
  return _test[_i][_qp] * val;
}

Real
DisloAnnihilation::computeQpJacobian()
{
  Real val;
  
  val = std::abs(_dislo_velocity[_qp][_slip_sys_index]);
  val *= _rho_annih[_qp];
  val *= _dc;
  val *= 4.0;
  
  return _test[_i][_qp] * _phi[_j][_qp] * val;
}

Real
DisloAnnihilation::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;	
	
  if (_rho_annih_coupled && jvar == _rho_annih_var)
  {
	  
	val = std::abs(_dislo_velocity[_qp][_slip_sys_index]);
	val *= _u[_qp];
    val *= _dc;
    val *= 4.0;
	
    return _test[_i][_qp] * _phi[_j][_qp] * val;

  }
  else {
	  
	return 0.0;
	
  } 
}
