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
  params.addCoupledVar("temp",303.0,"Temperature");  
  params.addParam<Real>("p_cs", 0.0,"Probability rate for cross slip in Christophe Depres thesis. "
                                    "This is the pre-factor that multiplies total screw "
									"density and gives rate of increase of curvature density. ");
  params.addParam<Real>("remain_rho_tol",0.000001,"Tolerance on small values of remain_rho_tot.");
  params.addParam<Real>("dislo_mobility",0.0,"Dislocation mobility");
  params.addParam<Real>("cross_slip_schmid_factor",0.0,"Ratio between Schmid factor of the "
                                                       "cross slip system and of the primary system");  
  params.addParam<Real>("tauIII",0.0,"Stage III resolved shear stress");
  params.addParam<Real>("dtauIII_dT",0.0,"Thermal coefficient of stage III resolved shear stress");
  params.addParam<Real>("reference_temperature",303.0,"reference temperature for stage III resolved shear stress");
  params.addParam<Real>("kB",0.0,"Boltzmann constant");
  params.addParam<Real>("Vact",0.0,"Activation volume");
  params.addParam<Real>("drho_cs_tol",0.0,"Upper limit of the cross slip rate of dislocation density");  
  params.addParam<Real>("R_cs",0.0,"Radius of cross slipped dislocation segment");
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
	_temp(coupledValue("temp")),
	_p_cs(getParam<Real>("p_cs")),
	_remain_rho_tol(getParam<Real>("remain_rho_tol")),
	_cssf(getParam<Real>("cross_slip_schmid_factor")),
	_tauIII(getParam<Real>("tauIII")),
    _dtauIII_dT(getParam<Real>("dtauIII_dT")),
    _reference_temperature(getParam<Real>("reference_temperature")),
	_kB(getParam<Real>("kB")),
	_Vact(getParam<Real>("Vact")),
    _drho_cs_tol(getParam<Real>("drho_cs_tol")),
    _R_cs(getParam<Real>("R_cs")),
	_tau_out(getMaterialProperty<std::vector<Real>>("tau_out")) // Resolved shear stress (signed)
{
}

Real
DoubleCrossSlip::computeQpResidual()
{
  Real val;
  Real remain_rho_tot;
  Real rss_cross_slip = 0.0; // resolved shear stress on the cross slip system
  Real temp = _temp[_qp]; // current temperature
  Real tauIII_T = 0.0;  // Stage III stress at temperature T
  Real drho_cs = 0.0; // cross slip rate of screw dislocation density

  remain_rho_tot = _rho_tot[_qp]*_rho_tot[_qp]
	             - _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp];
				 
  remain_rho_tot = std::max(remain_rho_tot,0.0);
  
  tauIII_T = _tauIII + _dtauIII_dT * (temp - _reference_temperature);
  
  // calculate resolved shear stress on the cross slip system
  // it is always positive here
  rss_cross_slip = _cssf * std::abs(_tau_out[_qp][_slip_sys_index]);
  
  // val is positive
  val = std::exp(((rss_cross_slip - tauIII_T) * _Vact) / (_kB * temp));
  val = _p_cs * val * std::sqrt(remain_rho_tot);

  drho_cs = val * _R_cs;

  if (drho_cs >= _drho_cs_tol) {

     val = _drho_cs_tol / _R_cs;

  }

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
  Real jac = 0.0;
  Real remain_rho_tot;
  Real rss_cross_slip = 0.0; // resolved shear stress on the cross slip system
  Real val;
  Real temp = _temp[_qp]; // current temperature
  Real tauIII_T = 0.0;  // Stage III stress at temperature T
  Real drho_cs = 0.0; // cross slip rate of screw dislocation density

  remain_rho_tot = _rho_tot[_qp]*_rho_tot[_qp]
	             - _rho_gnd_edge[_qp]*_rho_gnd_edge[_qp];
				 
  remain_rho_tot = std::max(remain_rho_tot,0.0);
				 
  tauIII_T = _tauIII + _dtauIII_dT * (temp - _reference_temperature);	
  
  // calculate resolved shear stress on the cross slip system
  // it is always positive here
  rss_cross_slip = _cssf * std::abs(_tau_out[_qp][_slip_sys_index]);
				 
  // val is positive
  val = std::exp(((rss_cross_slip - tauIII_T) * _Vact) / (_kB * temp));
  val = _p_cs * val * std::sqrt(remain_rho_tot);

  drho_cs = val * _R_cs;
  
  if (drho_cs >= _drho_cs_tol) { 
  
    // return 0 if rate is above limit because
	// rate is constant and does not depend on _rho_tot[_qp]
	// or on _rho_gnd_edge[_qp]

	return jac;

  }

  if (remain_rho_tot > _remain_rho_tol) { // check that denominator is not close to zero
	  
    if (_rho_tot_coupled && jvar == _rho_tot_var) {

      jac = - val
	      * (_rho_tot[_qp] / remain_rho_tot)
	      * _phi[_j][_qp] * _test[_i][_qp];
	
    } else if (_rho_gnd_edge_coupled && jvar == _rho_gnd_edge_var) {

      jac = val
          * (_rho_gnd_edge[_qp] / remain_rho_tot) 
	      * _phi[_j][_qp] * _test[_i][_qp];
	 
    }	  
	
  } // end check that denominator is not close to zero

  return jac; 
}
