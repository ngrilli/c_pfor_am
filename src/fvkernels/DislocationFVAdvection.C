// Nicolò Grilli
// Università di Bristol
// 22 Novembre 2025

#include "DislocationFVAdvection.h"
#include "Steady.h"
#include "FEProblemBase.h"

registerADMooseObject("MooseApp", DislocationFVAdvection);

InputParameters
DislocationFVAdvection::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params.addClassDescription(
      "Residual contribution from advection operator for finite volume method. "
      "This is a modification of FVAdvection "
      "where the velocity vector is imported from a material class "
      "that calculates dislocation velocity "
      "based on crystallographic directions and stress. ");
  params.addRequiredParam<Real>("velocity", "Constant velocity magnitude");
  params += Moose::FV::advectedInterpolationParameter();

  // We add the relationship manager here, this will select the right number of
  // ghosting layers depending on the chosen interpolation method
  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC | Moose::RelationshipManagerType::ALGEBRAIC |
          Moose::RelationshipManagerType::COUPLING,
      [](const InputParameters & obj_params, InputParameters & rm_params)
      { FVRelationshipManagerInterface::setRMParamsAdvection(obj_params, rm_params, 2); });

  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  return params;
}

DislocationFVAdvection::DislocationFVAdvection(const InputParameters & params)
  : FVFluxKernel(params), 
    _velocity(getParam<Real>("velocity")),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
    _screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
    _slip_sys_index(getParam<int>("slip_sys_index"))
{
  const bool need_more_ghosting =
      Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");
  if (need_more_ghosting && _tid == 0)
    // If we need more ghosting, then we are a second-order nonlinear limiting scheme whose stencil
    // is liable to change upon wind-direction change. Consequently we need to tell our problem that
    // it's ok to have new nonzeros which may crop-up after PETSc has shrunk the matrix memory
    getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")
        ->setErrorOnJacobianNonzeroReallocation(false);

  if (dynamic_cast<Steady *>(_app.getExecutioner()))
  {
    const MooseEnum not_available_with_steady("sou min_mod vanLeer quick venkatakrishnan");
    const std::string chosen_scheme =
        static_cast<std::string>(getParam<MooseEnum>("advected_interp_method"));
    if (not_available_with_steady.find(chosen_scheme) != not_available_with_steady.items().end())
      paramError("advected_interp_method",
                 "The given advected interpolation cannot be used with steady-state runs!");
  }
}

ADReal
DislocationFVAdvection::computeQpResidual()
{
  const auto state = determineState();
  const auto & limiter_time = _subproblem.isTransient()
                                  ? Moose::StateArg(1, Moose::SolutionIterationType::Time)
                                  : Moose::StateArg(1, Moose::SolutionIterationType::Nonlinear);

  RealVectorValue edge_velocity;
  RealVectorValue screw_velocity;
  
  // Allocate dislocation velocities based on slip systems index and dislocation character
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j) 
  {
    edge_velocity(j) = _edge_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // edge direction	  	
    screw_velocity(j) = - _screw_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // screw direction
  }

  const bool elem_is_upwind = _velocity * edge_velocity * _normal >= 0; // TO DO: edge vs screw
  const auto face = makeFace(*_face_info,
                             Moose::FV::limiterType(_advected_interp_method),
                             elem_is_upwind,
                             false,
                             &limiter_time);
  ADReal u_interface = _var(face, state);

  return _velocity * _normal * edge_velocity * u_interface;
}
