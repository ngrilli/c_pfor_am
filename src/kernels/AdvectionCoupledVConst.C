// Nicolo Grilli
// University of Bristol
// 6 Luglio 2021

#include "AdvectionCoupledVConst.h"
#include "SystemBase.h"

registerMooseObject("MooseApp", AdvectionCoupledVConst);

InputParameters
AdvectionCoupledVConst::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Coupled Advection. "
							 "Constant velocity. ");
  params.addRequiredCoupledVar("velocity", "Velocity vector");
  params.addCoupledVar("rho_coupled", 0.0, "Coupled dislocation density in the flux term.");
  return params;
}

AdvectionCoupledVConst::AdvectionCoupledVConst(const InputParameters & parameters)
  : Kernel(parameters),
    _velocity(coupledVectorValue("velocity")),
    _rho_coupled(coupledValue("rho_coupled")), // Coupled dislocation density in the flux term
    _rho_coupled_coupled(isCoupled("rho_coupled")),
    _rho_coupled_var(_rho_coupled_coupled ? coupled("rho_coupled") : 0)
{
}

Real
AdvectionCoupledVConst::negSpeedQp() const
{
  return -_grad_test[_i][_qp] * _velocity[_qp];
}

Real
AdvectionCoupledVConst::computeQpResidual()
{
  // This is the no-upwinded version
  // It gets called via Kernel::computeResidual()
  return negSpeedQp() * _rho_coupled[_qp];
}

Real
AdvectionCoupledVConst::computeQpJacobian()
{
  // This is the no-upwinded version
  // It gets called via Kernel::computeJacobian()
  return 0.0;
}

Real
AdvectionCoupledVConst::computeQpOffDiagJacobian(unsigned int jvar)
{
  // It gets called via Kernel::computeOffDiagJacobian()
  
  if (_rho_coupled_coupled && jvar == _rho_coupled_var)
  {
		
    return negSpeedQp() * _phi[_j][_qp];
	
  }
  else {
	  
	return 0.0;
	
  }
}

void
AdvectionCoupledVConst::computeResidual()
{
  Kernel::computeResidual();
}

void
AdvectionCoupledVConst::computeJacobian()
{
  Kernel::computeJacobian();
}

void
AdvectionCoupledVConst::computeOffDiagJacobian(const unsigned int jvar_num)
{
  Kernel::computeOffDiagJacobian(jvar_num);
}
