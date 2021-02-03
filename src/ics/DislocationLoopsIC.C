// Nicolò Grilli
// 1 Febbraio 2021
// National University of Singapore

#include "DislocationLoopsIC.h"

#include <cmath>

registerMooseObject("MooseApp", DislocationLoopsIC);

defineLegacyParams(DislocationLoopsIC);

InputParameters
DislocationLoopsIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addParam<std::vector<Real>>("centrex", "X coordinates of the centres of the loops");
  params.addParam<std::vector<Real>>("centrey", "Y coordinates of the centres of the loops");
  params.addParam<std::vector<Real>>("radii", "Radii of the loops");
  params.addParam<std::vector<Real>>("width", "Width of the loops");
  params.addParam<std::vector<Real>>("rho_max", "Max dislocation density in one loop");
  MooseEnum var_options("rhotot rhoedgegnd rhoscrewgnd qtot", "rhotot");
  params.addParam<MooseEnum>("variable_type",
                             var_options,
                             "Type of variable for the dislocation loop.");  
  params.addClassDescription("Initialize dislocation loops.");
  return params;
}

DislocationLoopsIC::DislocationLoopsIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _centrex(getParam<std::vector<Real>>("centrex")),
    _centrey(getParam<std::vector<Real>>("centrey")),
	_radii(getParam<std::vector<Real>>("radii")),
	_width(getParam<std::vector<Real>>("width")),
	_rho_max(getParam<std::vector<Real>>("rho_max")),
    _variable_type(getParam<MooseEnum>("variable_type"))	
{
}

// This function assumes x and y slip plane
// must be generalised to an arbitrary orientation of the slip plane
Real
DislocationLoopsIC::value(const Point & p)
{
	
  Real rhotot = 0.0;
  Real rhototloop; // temporary variable for dislocation loops
  Real rhoedgegnd = 0.0;
  Real rhoscrewgnd = 0.0;
  Real qtot = 0.0;
  Real R; // distance from centre, temporary variable
  Real radialvx; // Temporary X component of the unit radial vector
  Real radialvy; // Temporary Y component of the unit radial vector
  Real val;
	
  if (_centrex.size() <= 0 || _centrey.size() <= 0 || _radii.size() <= 0 || _width.size() <= 0)
    mooseError("Error in reading dislocation loops properties");

  for (unsigned int i = 0; i < _centrex.size(); ++i) {
	  
	R = std::sqrt((p(0)-_centrex[i])*(p(0)-_centrex[i])+(p(1)-_centrey[i])*(p(1)-_centrey[i]));
	
	rhototloop = std::exp(-(R-_radii[i])*(R-_radii[i])/(2.0*_width[i]*_width[i]));
	rhototloop *= (1.0/(2.506627 * _width[i]));
	rhototloop *= _rho_max[i];
	
	rhotot += rhototloop; 
	qtot = qtot + (rhototloop / _radii[i]);
	
	// Multiply total dislocation density by the unit radial vector
	// to obtain dislocation density vector
	if (R > 0.0001) {
	  radialvx = (p(0)-_centrex[i]) / R;
	  radialvy = (p(1)-_centrey[i]) / R;		
	} else {
	  radialvx = 0.0;
	  radialvy = 0.0;
	}

	rhoedgegnd = rhoedgegnd + rhototloop * radialvx;
	rhoscrewgnd = rhoscrewgnd + rhototloop * radialvy;
  }

  switch (_variable_type)
  {
    case 0: // rhotot
      val = rhotot;
      break;
    case 1: // rhoedgegnd
      val = rhoedgegnd;
	  break;
	case 2: // rhoscrewgnd
	  val = rhoscrewgnd;
	  break;
    case 3: // qtot
	  val = qtot;
      break;	
  }
  
  return val;

}
