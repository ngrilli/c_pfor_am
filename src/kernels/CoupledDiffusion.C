// Nicolò Grilli
// National University of Singapore
// 2 Febbraio 2021

#include "CoupledDiffusion.h"

registerMooseObject("MooseApp", CoupledDiffusion);

InputParameters
CoupledDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Residual of: coef * nabla^2 (coupled_variable)");
  params.addRequiredCoupledVar("coupled_variable", "Coupled variable for Laplacian operator");
  params.addCustomTypeParam("coef", 0.0, "CoefficientType", "The coefficient of diffusion");
  params.declareControllable("coef");
  return params;
}

CoupledDiffusion::CoupledDiffusion(const InputParameters & parameters)
  : Kernel(parameters),
    _velocity_vector(coupledGradient("coupled_variable")),
    _is_var_coupled(isCoupled("coupled_variable")),
    _coupled_var(_is_var_coupled ? coupled("coupled_variable") : 0),
	_coef(getParam<Real>("coef"))
{
}

Real
CoupledDiffusion::computeQpResidual()
{
  return _coef * (_grad_test[_i][_qp] * _velocity_vector[_qp]);
}

Real
CoupledDiffusion::computeQpJacobian()
{
  // Residual depends only on the coupled variable
  return 0.0;
}

Real
CoupledDiffusion::computeQpOffDiagJacobian(unsigned int jvar)
{
	
  if (_is_var_coupled && jvar == _coupled_var)
  {
	
	return _coef * (_grad_test[_i][_qp] * _grad_phi[_j][_qp]);

  }
  else {
	  
	return 0.0;
	
  } 	
}

