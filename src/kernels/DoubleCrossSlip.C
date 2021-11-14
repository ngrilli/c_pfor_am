// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// University of Bristol
// 13 Novembre 2021

#include "DoubleCrossSlip.h"

registerMooseObject("MooseApp", DoubleCrossSlip);

template <>
InputParameters
validParams<DoubleCrossSlip>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Dislocation double cross slip. "
                             "This kernel is meant to be applied to the curvature equation. "
							 "It represents the creation of new curvature after a "
							 "double slip process and bow out of the dislocation loop "
							 "on a parallel slip plane. ");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip velocity, "
							   "for instance from 0 to 11 for FCC."); 
  params.addCoupledVar("rho_gnd_edge", 0.0, "Edge GND dislocation density");
  params.addCoupledVar("rho_gnd_screw", 0.0, "Screw GND dislocation density");
  params.addCoupledVar("rho_tot", 0.0, "Total dislocation density: rho_t.");               
  params.addParam<Real>("p_cs", 0.0,"Probability rate for cross slip in El-Azab 2016 paper. "
                                    "This is the pre-factor that multiplies total screw "
									"density and gives rate of increase of curvature density. ");
  params.addParam<Real>("remain_rho_tol",0.000001,"Tolerance on small values of remain_rho_tot.");									
  return params;
}

DoubleCrossSlip::DoubleCrossSlip(const InputParameters & parameters)
  : Kernel(parameters),
	_slip_sys_index(getParam<int>("slip_sys_index")),
    _rho_gnd_edge(coupledValue("rho_gnd_edge")), // GND dislocation density: rho_x or rho_y for edge or screw
    _rho_gnd_edge_coupled(isCoupled("rho_gnd_edge")),
    _rho_gnd_edge_var(_rho_gnd_edge_coupled ? coupled("rho_gnd_edge") : 0),
    _rho_gnd_screw(coupledValue("rho_gnd_screw")), // GND dislocation density: rho_x or rho_y for edge or screw
    _rho_gnd_screw_coupled(isCoupled("rho_gnd_screw")),
    _rho_gnd_screw_var(_rho_gnd_screw_coupled ? coupled("rho_gnd_screw") : 0),	
    _rho_tot(coupledValue("rho_tot")), // Total dislocation density: rho_t
    _rho_tot_coupled(isCoupled("rho_tot")),
    _rho_tot_var(_rho_tot_coupled ? coupled("rho_tot") : 0),
	_p_cs(getParam<Real>("p_cs")),
	_remain_rho_tol(getParam<Real>("remain_rho_tol")),
	_dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")) // Velocity value (signed)
{
}

Real
DoubleCrossSlip::computeQpResidual()
{
  Real val;
  Real remain_rho_tot;

  remain_rho_tot = _rho_tot[_qp]*_rho_tot[_qp]
	             - _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp];
				 
  remain_rho_tot = std::max(remain_rho_tot,0.0);

  // val is positive
  val = _p_cs * std::abs(_dislo_velocity[_qp][_slip_sys_index]) * std::sqrt(remain_rho_tot);

  return - _test[_i][_qp] * val; // minus sign because this is in the LHS of the equation
}

// Residual does not depend on the curvature density
Real
DoubleCrossSlip::computeQpJacobian()
{
  return 0.0;
}

Real
DoubleCrossSlip::computeQpOffDiagJacobian(unsigned int jvar)
{	  
  Real jac = 0;
  Real remain_rho_tot;

  remain_rho_tot = _rho_tot[_qp]*_rho_tot[_qp]
	             - _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp];

  if (remain_rho_tot > _remain_rho_tol) { // check that denominator is not close to zero
	  
    if (_rho_tot_coupled && jvar == _rho_tot_var) {

      jac = - _p_cs * std::abs(_dislo_velocity[_qp][_slip_sys_index]) 
	      * (_rho_tot[_qp] / std::sqrt(remain_rho_tot))
	      * _phi[_j][_qp] * _test[_i][_qp];
	
    } else if (_rho_gnd_edge_coupled && jvar == _rho_gnd_edge_var) {

      jac = _p_cs * std::abs(_dislo_velocity[_qp][_slip_sys_index])
          * (_rho_gnd_edge[_qp] / std::sqrt(remain_rho_tot)) 
	      * _phi[_j][_qp] * _test[_i][_qp];
	 
    }	  
	
  } // end check that denominator is not close to zero

  return jac; 
}
