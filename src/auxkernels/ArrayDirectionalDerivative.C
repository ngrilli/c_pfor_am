// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// University of Bristol
// 31 Marzo 2022

#include "ArrayDirectionalDerivative.h"

registerMooseObject("MooseApp", ArrayDirectionalDerivative);

InputParameters
ArrayDirectionalDerivative::validParams()
{
  InputParameters params = ArrayAuxKernel::validParams();
  params.addClassDescription(
      "Calculate directional derivative along edge and screw dislocation propagation direction."
	  "This AuxKernel applies to a vector auxiliary variable");
  params.addRequiredCoupledVar("gradient_variable",
                               "The vector variable from which to compute the directional derivative");
  MooseEnum dislo_character("edge screw", "edge");
  params.addRequiredParam<MooseEnum>("dislo_character",
                                     dislo_character,
                                     "Character of dislocations: edge or screw.");
  return params;
}

ArrayDirectionalDerivative::ArrayDirectionalDerivative(const InputParameters & parameters)
  : ArrayAuxKernel(parameters),
    _grad_variable(coupledArrayGradient("gradient_variable")),
    _edge_slip_direction(getMaterialProperty<std::vector<Real>>("edge_slip_direction")), // Edge velocity direction
	_screw_slip_direction(getMaterialProperty<std::vector<Real>>("screw_slip_direction")), // Screw velocity direction
	_dislo_character(getParam<MooseEnum>("dislo_character").getEnum<DisloCharacter>())
{
}

// Directional derivative along the edge or screw slip direction
// It must be calculated for all slip systems
RealEigenVector
ArrayDirectionalDerivative::computeValue()
{
  RealEigenVector val(_var.count());

  // Temporary vector to store the gradient components
  // of the variable on a specific slip system
  //RealVectorValue temp_grad_variable;

  _velocity.resize(3, 0.0);
  
  // Find dislocation velocity based on slip systems index and dislocation character
  // _velocity here is a temporary vector that change when considering each slip system
  for (unsigned int i = 0; i < _var.count(); ++i){
	  
  // check if edge or screw direction for slip system with index i
  switch (_dislo_character)
  {
    case DisloCharacter::edge:
	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	  {
	    _velocity[j] = _edge_slip_direction[_qp][i * LIBMESH_DIM + j]; // edge direction	  
	  }
	  break;
	case DisloCharacter::screw:
	  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	  {
		// note that the definition of _screw_slip_direction in FiniteStrainCrystalPlasticityDislo
		// and CrystalPlasticityDislocationUpdate
		// is -y, because +x is _edge_slip_direction and +z is slip plane normal
		// but derivative must be taken along +y
		// therefore a sign change is needed
	    _velocity[j] = - _screw_slip_direction[_qp][i * LIBMESH_DIM + j]; // screw direction	  
	  }	
	  break;
  }  

  // Calculate output for slip system i
  val(i) = 0.0;
  
  // _grad_variable[_qp] is an array with number of component
  // given by _var.count() * LIBMESH_DIM
  // the order of the components is, for instance, in 3D
  // the first _var.count() variables are the derivatives along x of
  // the components of _var
  // the next _var.count() variables are the derivatives along y of
  // the components of _var
  // the next _var.count() variables are the derivatives along z of
  // the components of _var
    
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j) {
    val(i) += _grad_variable[_qp](j*_var.count()+i) * _velocity[j];
  }
	
  } // end of iteration over slip systems
  
  return val;
}
