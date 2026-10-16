// Nicolò Grilli
// Università di Bristol
// 14 Ottobre 2026

#include "CurlCurlComponent.h"

registerMooseObject("c_pfor_amApp", CurlCurlComponent);

InputParameters
CurlCurlComponent::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Set the kernel variable to a specified component of the curl curl of a vector made of " 
                             "three variables Hx, Hy, Hz, where the two variables involved in the calculations "
                             "are considered coupled variables, while the variable not involved is considered "
                             "an independent variable. The weak form of the curl curl is calculated as follows: "
                             "alpha * curl H * curl test");
  params.addRequiredCoupledVar("H1", "Coupled variable whose derivative has positive sign in the curl component");
  params.addRequiredCoupledVar("H2", "Coupled variable whose derivative has negative sign in the curl component");
  params.addRequiredParam<unsigned int>("component","Component of the curl curl of Hx, Hy, Hz calculated by this kernel");
  params.addParam<MaterialPropertyName>("resistivity", "resistivity", "Name of the material property for current-dependent resistivity");
  return params;
}

CurlCurlComponent::CurlCurlComponent(const InputParameters & parameters)
  : Kernel(parameters),
    _H1_var(coupled("H1")),
    _H2_var(coupled("H2")),
    _grad_H1(coupledGradient("H1")),
    _grad_H2(coupledGradient("H2")),
    _component(getParam<unsigned int>("component")),
    _resistivity(getMaterialProperty<Real>("resistivity"))
{
  if (_component >= LIBMESH_DIM)
    paramError("component", "Component too large for LIBMESH_DIM");
}

// Curl of H: output is the same vector in all cases but different cases
// have different variable names meaning different components of H
void
CurlCurlComponent::computeCurlH(Real & CurlHx, Real & CurlHy, Real & CurlHz)
{
  // Cases:
  // component 0: u=Hx, H1=Hz, H2=Hy
  // component 1: u=Hy, H1=Hx, H2=Hz
  // component 2: u=Hz, H1=Hy, H2=Hx

  if (_component == 0)
  {
    CurlHx = _grad_H1[_qp](1) - _grad_H2[_qp](2);
    CurlHy = _grad_u[_qp](2) - _grad_H1[_qp](0);
    CurlHz = _grad_H2[_qp](0) - _grad_u[_qp](1);
  }
  else if (_component == 1)
  {
    CurlHx = _grad_H2[_qp](1) - _grad_u[_qp](2);
    CurlHy = _grad_H1[_qp](2) - _grad_H2[_qp](0);
    CurlHz = _grad_u[_qp](0) - _grad_H1[_qp](1);
  }
  else if (_component == 2)
  {
    CurlHx = _grad_u[_qp](1) - _grad_H1[_qp](2);
    CurlHy = _grad_H2[_qp](2) - _grad_u[_qp](0);
    CurlHz = _grad_H1[_qp](0) - _grad_H2[_qp](1);
  }
}

Real
CurlCurlComponent::computeQpResidual()
{
  Real CurlHx, CurlHy, CurlHz;
  computeCurlH(CurlHx, CurlHy, CurlHz);

  if (_component == 0) // component 0: u=Hx, H1=Hz, H2=Hy
    return _resistivity[_qp] * (CurlHy * _grad_test[_i][_qp](2) - CurlHz * _grad_test[_i][_qp](1));
  else if (_component == 1) // component 1: u=Hy, H1=Hx, H2=Hz
    return _resistivity[_qp] * (CurlHz * _grad_test[_i][_qp](0) - CurlHx * _grad_test[_i][_qp](2));
  else if (_component == 2) // component 2: u=Hz, H1=Hy, H2=Hx
    return _resistivity[_qp] * (CurlHx * _grad_test[_i][_qp](1) - CurlHy * _grad_test[_i][_qp](0));
  else
    return 0.0;
}

Real
CurlCurlComponent::computeQpJacobian()
{
  if (_component == 0) // component 0: u=Hx, H1=Hz, H2=Hy
    return _resistivity[_qp] * (_grad_phi[_j][_qp](2) * _grad_test[_i][_qp](2) + _grad_phi[_j][_qp](1) * _grad_test[_i][_qp](1));
  else if (_component == 1) // component 1: u=Hy, H1=Hx, H2=Hz
    return _resistivity[_qp] * (_grad_phi[_j][_qp](0) * _grad_test[_i][_qp](0) + _grad_phi[_j][_qp](2) * _grad_test[_i][_qp](2));
  else if (_component == 2) // component 2: u=Hz, H1=Hy, H2=Hx
    return _resistivity[_qp] * (_grad_phi[_j][_qp](1) * _grad_test[_i][_qp](1) + _grad_phi[_j][_qp](0) * _grad_test[_i][_qp](0));
  else
    return 0.0;
}

Real
CurlCurlComponent::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0) // component 0: u=Hx, H1=Hz, H2=Hy
  {
    if (jvar == _H1_var) // Take derivative of CurlHy
      return -_resistivity[_qp] * _grad_phi[_j][_qp](0) * _grad_test[_i][_qp](2);
    else if (jvar == _H2_var) // Take derivative of CurlHz
      return -_resistivity[_qp] * _grad_phi[_j][_qp](0) * _grad_test[_i][_qp](1);
  }
  else if (_component == 1) // component 1: u=Hy, H1=Hx, H2=Hz
  {
    if (jvar == _H1_var) // Take derivative of CurlHz
      return -_resistivity[_qp] * _grad_phi[_j][_qp](1) * _grad_test[_i][_qp](0);
    else if (jvar == _H2_var) // Take derivative of CurlHx
      return -_resistivity[_qp] * _grad_phi[_j][_qp](1) * _grad_test[_i][_qp](2);
  }
  else if (_component == 2) // component 2: u=Hz, H1=Hy, H2=Hx
  {
    if (jvar == _H1_var) // Take derivative of CurlHx
      return -_resistivity[_qp] * _grad_phi[_j][_qp](2) * _grad_test[_i][_qp](1);
    else if (jvar == _H2_var) // Take derivative of CurlHy
      return -_resistivity[_qp] * _grad_phi[_j][_qp](2) * _grad_test[_i][_qp](0);
  }
  return 0.0;
}
