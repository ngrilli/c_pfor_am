// Nicolo Grilli
// National University of Singapore
// 26 Febbraio 2021

#include "DirectionalDerivative.h"

registerMooseObject("MooseApp", DirectionalDerivative);

defineLegacyParams(DirectionalDerivative);

InputParameters
DirectionalDerivative::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Calculate directional derivative along edge and screw dislocation propagation direction.");
  params.addRequiredCoupledVar("gradient_variable",
                               "The variable from which to compute the directional derivative");
  params.addRequiredParam<int>("slip_sys_index", "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  return params;
}

DirectionalDerivative::DirectionalDerivative(const InputParameters & parameters)
  : AuxKernel(parameters),
    _gradient(coupledGradient("gradient_variable")),
    _slip_sys_index(getParam<int>("slip_sys_index")),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>())
{
}

// Directional derivative along the edge or screw slip direction
Real
DirectionalDerivative::computeValue()
{
  Real val = 0.0;	
	
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
  
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j) {
	val += _gradient[_qp](j) * _velocity[j];  
  }
	
  return val;
}
