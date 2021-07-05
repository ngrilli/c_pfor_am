// Nicolò Grilli
// University of Bristol
// 5 Luglio 2021

#include "FVCoupledAdvection.h"

registerADMooseObject("MooseApp", FVCoupledAdvection);

InputParameters
FVCoupledAdvection::validParams()
{
  InputParameters params = FVCoupledFluxKernel::validParams();
  params.addClassDescription(
      "Residual contribution from advection operator for finite volume method "
	  "of a coupled variable.");
  params.addRequiredParam<RealVectorValue>("velocity", "Constant advection velocity");
  MooseEnum advected_interp_method("average upwind", "upwind");

  params.addParam<MooseEnum>("advected_interp_method",
                             advected_interp_method,
                             "The interpolation to use for the advected quantity. Options are "
                             "'upwind' and 'average', with the default being 'upwind'.");
  return params;
}

FVCoupledAdvection::FVCoupledAdvection(const InputParameters & params)
  : FVCoupledFluxKernel(params), _velocity(getParam<RealVectorValue>("velocity"))
{
  using namespace Moose::FV;

  const auto & advected_interp_method = getParam<MooseEnum>("advected_interp_method");
  if (advected_interp_method == "average")
    _advected_interp_method = InterpMethod::Average;
  else if (advected_interp_method == "upwind")
    _advected_interp_method = InterpMethod::Upwind;
  else
    mooseError("Unrecognized interpolation type ",
               static_cast<std::string>(advected_interp_method));
}

ADReal
FVCoupledAdvection::computeQpResidual()
{
  ADReal u_interface;
  interpolate(
      _advected_interp_method, u_interface, _u_elem[_qp], _u_neighbor[_qp], _velocity, *_face_info);
  return _normal * _velocity * u_interface;
}
