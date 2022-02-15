// Nicolò Grilli
// University of Bristol
// 15 Febbraio 2022

#include "ACInterfaceSlipPlaneFracture.h"

registerMooseObject("PhaseFieldApp", ACInterfaceSlipPlaneFracture);

InputParameters
ACInterfaceSlipPlaneFracture::validParams()
{
  InputParameters params = ACInterface::validParams();
  params.addClassDescription("Gradient energy Allen-Cahn Kernel where crack propagation along weak"
                             "cleavage plane is preferred."
							 "Considers cleavage plane anisotropy in the crack propagation "
							 "the cleavage plane is based on a pre-selected slip plane.");
  params.addRequiredParam<Real>(
      "beta_penalty",
      "penalty to penalize fracture on planes not normal to one cleavage plane normal which is "
      "normal to weak cleavage plane. Setting beta=0 results in isotropic damage.");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip normal "
							   "used for the cleavage plane.");
  return params;
}

ACInterfaceSlipPlaneFracture::ACInterfaceSlipPlaneFracture(const InputParameters & parameters)
  : ACInterface(parameters),
    _beta_penalty(getParam<Real>("beta_penalty")),
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_slip_plane_normals(getMaterialProperty<std::vector<Real>>("slip_plane_normals"))
{
}

Real
ACInterfaceSlipPlaneFracture::betaNablaPsi()
{
  // assign cleavage plane based on slip system	
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j) {
    _cleavage_plane_normal(j) = _slip_plane_normals[_qp][_slip_sys_index * LIBMESH_DIM + j];
  }	
	
  return _beta_penalty * _L[_qp] * _kappa[_qp] * (_grad_u[_qp] * _cleavage_plane_normal) *
         (_grad_test[_i][_qp] * _cleavage_plane_normal);
}

Real
ACInterfaceSlipPlaneFracture::computeQpResidual()
{
  return (1 + _beta_penalty) * _grad_u[_qp] * kappaNablaLPsi() - betaNablaPsi();
}

Real
ACInterfaceSlipPlaneFracture::computeQpJacobian()
{
  // assign cleavage plane based on slip system	
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j) {
    _cleavage_plane_normal(j) = _slip_plane_normals[_qp][_slip_sys_index * LIBMESH_DIM + j];;
  }		
	
  /// dsum is the derivative \f$ \frac\partial{\partial \eta} \left( \nabla
  /// (L\psi) \right) \f$
  RealGradient dsum =
      (_dkappadop[_qp] * _L[_qp] + _kappa[_qp] * _dLdop[_qp]) * _phi[_j][_qp] * _grad_test[_i][_qp];

  /// compute the derivative of the gradient of the mobility
  if (_variable_L)
  {
    RealGradient dgradL =
        _grad_phi[_j][_qp] * _dLdop[_qp] + _grad_u[_qp] * _phi[_j][_qp] * _d2Ldop2[_qp];

    for (unsigned int i = 0; i < _n_args; ++i)
      dgradL += (*_gradarg[i])[_qp] * _phi[_j][_qp] * (*_d2Ldargdop[i])[_qp];

    dsum += (_kappa[_qp] * dgradL + _dkappadop[_qp] * _phi[_j][_qp] * gradL()) * _test[_i][_qp];
  }

  return (1 + _beta_penalty) * _grad_phi[_j][_qp] * kappaNablaLPsi() + _grad_u[_qp] * dsum -
         _beta_penalty * _L[_qp] * _kappa[_qp] * (_grad_u[_qp] * _cleavage_plane_normal) *
             (_grad_phi[_j][_qp] * _cleavage_plane_normal);
}
