// Nicolo Grilli
// National University of Singapore
// 18 Gennaio 2021

#include "ConservativeAdvectionCoupled.h"
#include "SystemBase.h"
#include "libmesh/utility.h"

registerMooseObject("MooseApp", ConservativeAdvectionCoupled);

defineLegacyParams(ConservativeAdvectionCoupled);

InputParameters
ConservativeAdvectionCoupled::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u_coupl$ which in its weak "
                             "form is given by: $(-\\nabla \\psi_i, \\vec{v} u_coupl)$. "
							 "Velocity \vec{v} is taken as material property."
							 "u_coupl is a coupled variable");
  params.addCoupledVar("rho_coupled", 0.0, "Coupled dislocation density in the flux term.");
  MooseEnum upwinding_type("none full", "none");
  params.addParam<MooseEnum>("upwinding_type",
                             upwinding_type,
                             "Type of upwinding used.  None: Typically results in overshoots and "
                             "undershoots, but numerical diffusion is minimized.  Full: Overshoots "
                             "and undershoots are avoided, but numerical diffusion is large");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  MooseEnum dislo_sign("positive negative", "positive");
  params.addRequiredParam<MooseEnum>("dislo_sign",
                                     dislo_sign,
                                     "Sign of dislocations.");
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  return params;
}

ConservativeAdvectionCoupled::ConservativeAdvectionCoupled(const InputParameters & parameters)
  : Kernel(parameters),
    _rho_coupled(coupledValue("rho_coupled")), // Coupled dislocation density in the flux term
    _rho_coupled_coupled(isCoupled("rho_coupled")),
    _rho_coupled_var(_rho_coupled_coupled ? coupled("rho_coupled") : 0),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
    _dislo_velocity(getMaterialProperty<std::vector<Real>>("dislo_velocity")), // Velocity value (signed)
    _upwinding(getParam<MooseEnum>("upwinding_type").getEnum<UpwindingType>()),
	_slip_sys_index(getParam<int>("slip_sys_index")),
	_dislo_sign(getParam<MooseEnum>("dislo_sign").getEnum<DisloSign>()),
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>()),
	_u_coupled_nodal(getVar("rho_coupled", 0)->dofValues()), // nodal values of coupled dislocation density
    //_u_nodal(_var.dofValues()), // is this useful?
    _upwind_node(0),
    _dtotal_mass_out(0)
{
}

Real
ConservativeAdvectionCoupled::negSpeedQp()
{
  Real edge_sign;

  switch (_dislo_sign)
  {
    case DisloSign::positive:
      edge_sign = 1.0;
      break;
    case DisloSign::negative:
      edge_sign = -1.0;
      break;
  }  
	
  _velocity.resize(3, 0.0);
  
  // Find dislocation velocity based on slip systems index and dislocation character
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	  {
	    _velocity[j] = _edge_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // edge direction	  
	  }
	  break;
	case DisloCharacter::screw:
	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	  {
		// note that the definition of _screw_slip_direction in FiniteStrainCrystalPlasticityDislo
		// is -y, because +x is _edge_slip_direction and +z is slip plane normal
		// but derivative must be taken along +y
		// therefore a sign change is needed
	    _velocity[j] = - _screw_slip_direction[_qp][_slip_sys_index * LIBMESH_DIM + j]; // screw direction	  
	  }	
	  break;
  }
	
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  {
	_velocity[j] *= _dislo_velocity[_qp][_slip_sys_index]; // velocity value (signed)
	_velocity[j] *= edge_sign; // positive or negative dislocation
  }
  
  return -_grad_test[_i][_qp] * RealVectorValue(_velocity[0],_velocity[1],_velocity[2]);
}

Real
ConservativeAdvectionCoupled::computeQpResidual()
{
  // This is the no-upwinded version
  // It gets called via Kernel::computeResidual()
  return negSpeedQp() * _rho_coupled[_qp];
}

Real
ConservativeAdvectionCoupled::computeQpJacobian()
{
  // This is the no-upwinded version
  // It gets called via Kernel::computeJacobian()
  return 0.0;
}

Real
ConservativeAdvectionCoupled::computeQpOffDiagJacobian(unsigned int jvar)
{
  // This is the no-upwinded version
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
ConservativeAdvectionCoupled::computeResidual()
{
  switch (_upwinding)
  {
    case UpwindingType::none:
      Kernel::computeResidual();
      break;
    case UpwindingType::full:
      fullUpwind(JacRes::CALCULATE_RESIDUAL);
      break;
  }
}

void
ConservativeAdvectionCoupled::computeJacobian()
{
  // Always return zero in this case
  Kernel::computeJacobian();
}

void
ConservativeAdvectionCoupled::computeOffDiagJacobian(const unsigned int jvar_num)
{
  switch (_upwinding)
  {
    case UpwindingType::none:
      Kernel::computeOffDiagJacobian(jvar_num);
      break;
    case UpwindingType::full:
      fullUpwind(JacRes::CALCULATE_JACOBIAN);
      break;
  }	
}

void
ConservativeAdvectionCoupled::fullUpwind(JacRes res_or_jac)
{
  // The number of nodes in the element
  const unsigned int num_nodes = _test.size();

  // Even if we are computing the Jacobian we still need to compute the outflow from each node to
  // see which nodes are upwind and which are downwind
  prepareVectorTag(_assembly, _var.number());

  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
    prepareMatrixTag(_assembly, _var.number(), getVar("rho_coupled", 0)->number());

  // Compute the outflux from each node and store in _local_re
  // If _local_re is positive at the node, mass (or whatever the Variable represents) is flowing out
  // of the node
  _upwind_node.resize(num_nodes);
  for (_i = 0; _i < num_nodes; ++_i)
  {
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * negSpeedQp();
    _upwind_node[_i] = (_local_re(_i) >= 0.0);
  }

  // Variables used to ensure mass conservation
  Real total_mass_out = 0.0;
  Real total_in = 0.0;
  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
    _dtotal_mass_out.assign(num_nodes, 0.0);

  for (unsigned int n = 0; n < num_nodes; ++n)
  {
    if (_upwind_node[n])
    {
      if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
      {
        if (_test.size() == _phi.size())
          /* u_coupl at node=n depends only on the u_coupl at node=n, by construction.  For
           * linear-lagrange variables, this means that Jacobian entries involving the derivative
           * will only be nonzero for derivatives wrt coupled variable at node=n.  Hence the
           * (n, n) in the line below.  The above "if" statement catches other variable types
           * (eg constant monomials)
           */
          _local_ke(n, n) += _local_re(n); // off-diagonal Jacobian

        _dtotal_mass_out[n] += _local_ke(n, n);
      }
      _local_re(n) *= _u_coupled_nodal[n];
      total_mass_out += _local_re(n);
    }
    else                        // downwind node
      total_in -= _local_re(n); // note the -= means the result is positive
  }

  // Conserve mass over all phases by proportioning the total_mass_out mass to the inflow nodes,
  // weighted by their local_re values
  for (unsigned int n = 0; n < num_nodes; ++n)
  {
    if (!_upwind_node[n]) // downwind node
    {
      if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
        for (_j = 0; _j < _phi.size(); _j++)
          _local_ke(n, _j) += _local_re(n) * _dtotal_mass_out[_j] / total_in;
      _local_re(n) *= total_mass_out / total_in;
    }
  }

  // Add the result to the residual and jacobian
  if (res_or_jac == JacRes::CALCULATE_RESIDUAL)
  {
    accumulateTaggedLocalResidual();

    if (_has_save_in)
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (const auto & var : _save_in)
        var->sys().solution().add_vector(_local_re, var->dofIndices());
    }
  }

  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
  {
    accumulateTaggedLocalMatrix();

    if (_has_diag_save_in)
    {
      unsigned int rows = _local_ke.m();
      DenseVector<Number> diag(rows);
      for (unsigned int i = 0; i < rows; i++)
        diag(i) = _local_ke(i, i);

      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (const auto & var : _diag_save_in)
        var->sys().solution().add_vector(diag, var->dofIndices());
    }
  }
}
