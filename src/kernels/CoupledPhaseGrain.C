// Fernando Valiente Dies
// ANSTO
// Nicolò Grilli
// Parsa Esmati
// Università di Bristol
// 17 Dicembre 2023

#include "CoupledPhaseGrain.h"

registerMooseObject("c_pfor_amApp", CoupledPhaseGrain);

InputParameters
CoupledPhaseGrain::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Interaction between phases and grain orientations for zeta term. ");
  params.addRequiredCoupledVar("v",
                               "Array of coupled order parameter names for all the order parameters");
  params.addParam<Real>("A",0.0,"Prefactor. ");
  return params;
}

CoupledPhaseGrain::CoupledPhaseGrain(const InputParameters & parameters)
  : Kernel(parameters), 
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_var(coupledIndices("v")),
    _A(getParam<Real>("A"))
{
}

// Assign values of the order parameters to a vector
std::vector<Real>
CoupledPhaseGrain::assignOps()
{
  std::vector<Real> all_ops(_op_num);
  for (unsigned int i = 0; i < _op_num; ++i)
    all_ops[i] = (*_vals[i])[_qp];

  return all_ops;
}

Real
CoupledPhaseGrain::computeQpResidual()
{
  std::vector<Real> all_ops(_op_num);
  all_ops = assignOps();
  
  // Sum the squares of all the order parameters
  Real SumOPj = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
    SumOPj += all_ops[i] * all_ops[i];

  return _test[_i][_qp] * _A * (1.0 - _u[_qp]) * SumOPj;
}

// Derivative with respect to zeta
Real
CoupledPhaseGrain::computeQpJacobian()
{
  std::vector<Real> all_ops(_op_num);
  all_ops = assignOps();
  
  // Sum the squares of all the order parameters
  Real SumOPj = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
    SumOPj += all_ops[i] * all_ops[i];	
	
  return (-1.0) * _test[_i][_qp] * _A * _phi[_j][_qp] * SumOPj;
}

Real
CoupledPhaseGrain::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real jac;
	
  for (unsigned int i = 0; i < _op_num; ++i) {
    if (jvar == _vals_var[i]) { // only the current j is selected in the sum
		
      jac = _test[_i][_qp] * _A * (1.0 - _u[_qp]) * _phi[_j][_qp] * 2.0 * (*_vals[i])[_qp];
		
	}	  
  }

  return jac;    	
}
