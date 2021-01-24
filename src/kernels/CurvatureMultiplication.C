// Nicolò Grilli
// National University of Singapore
// 24 Gennaio 2021

// Last term on the right hand side of equation (1) in
// Stefan Sandfeld and Michael Zaiser
// Pattern formation in a minimal model of continuum
// dislocation plasticity
// Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp)

#include "CurvatureMultiplication.h"

registerMooseObject("MooseApp", CurvatureMultiplication);

template <>
InputParameters
validParams<CurvatureMultiplication>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Total dislocation multiplication rate "
                             "Due to total dislocation curvature.");
  params.addCoupledVar("curvature", 0.0, "Total curvature density q_t.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip velocity, "
							   "for instance from 0 to 11 for FCC.");							   
  return params;
}

CurvatureMultiplication::CurvatureMultiplication(const InputParameters & parameters)
  : Kernel(parameters),
    _curvature(coupledValue("curvature")),
	_slip_sys_index(getParam<int>("slip_sys_index")),
    _curvature_coupled(isCoupled("curvature")),
    _curvature_var(_curvature_coupled ? coupled("curvature") : 0),
	_dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")) // Velocity value (signed)
{
}

Real
CurvatureMultiplication::computeQpResidual()
{
  Real val;
  
  val = _dislo_velocity[_qp][_slip_sys_index];
  val *= _curvature[_qp];
  val *= -1.0;
  
  return _test[_i][_qp] * val;
}

Real
CurvatureMultiplication::computeQpJacobian()
{
  return 0.0;
}

Real
CurvatureMultiplication::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;	
	
  if (_curvature_coupled && jvar == _curvature_var)
  {
	  
    val = _dislo_velocity[_qp][_slip_sys_index];
    val *= -1.0;
	
    return _test[_i][_qp] * _phi[_j][_qp] * val;

  }
  else {
	  
	return 0.0;
	
  } 
}
