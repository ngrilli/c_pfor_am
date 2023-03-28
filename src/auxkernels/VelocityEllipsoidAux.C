// Fernando Valiente Dies
// ANSTO
// Nicolò Grilli
// University of Bristol
// 22 Marzo 2023

#include "VelocityEllipsoidAux.h"

#include "Function.h"

registerMooseObject("MooseApp", VelocityEllipsoidAux);

InputParameters
VelocityEllipsoidAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("This AuxKernel increases an AuxVariable from 0 to 1 if the qp is inside "
                                             "an ellipsoid that is moving according to paths defined by velocity and initial coordinates. "
											 "It can be applied to a level set variable "
                                             "to simulate the material deposition during wire arc additive manufacturing (WAAM) "
											 "together with ActDeactElementsCoupled. "
											 "It uses a postprocessor to start the next path to simulate a temperature controlled WAAM deposition process. ");
  params.addParam<Real>("low_level_set_var", 0.0, "The lowest value of the level set variable.");
  params.addParam<Real>("high_level_set_var", 1.0, "The highest value of the level set variable.");
  params.addRequiredParam<Real>("rx", "effective transverse ellipsoid radius");
  params.addRequiredParam<Real>("ry", "effective longitudinal ellipsoid radius");
  params.addRequiredParam<Real>("rz", "effective depth ellipsoid radius");
  params.addRequiredParam<RealVectorValue>("velocity", "Velocity vector of the material deposition source. ");
  
  // Every time the postprocessor condition is satisfied, the material deposition source is moved to the next set of coordinates
  params.addRequiredParam<std::vector<Real>>("init_x_coords", "Initial values of x coordinates of the material deposition source");
  params.addRequiredParam<std::vector<Real>>("init_y_coords", "Initial values of y coordinates of the material deposition source");
  params.addRequiredParam<std::vector<Real>>("init_z_coords", "Initial values of z coordinates of the material deposition source");
  
  params.addRequiredParam<std::vector<Real>>("scan_length","Total length during one scan. "
                                                           "After this length the material deposition is switched off. ");
  
  
  params.addParam<Real>("level_set_activation_threshold", 0.5, "Threshold value of the ellipsoid function "
                                                                                                 "that activates the level set.");	  
  return params;
}

VelocityEllipsoidAux::VelocityEllipsoidAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _low_level_set_var(getParam<Real>("low_level_set_var")),
	_high_level_set_var(getParam<Real>("high_level_set_var")),
    _rx(getParam<Real>("rx")),
    _ry(getParam<Real>("ry")),
    _rz(getParam<Real>("rz")),
    _velocity(getParam<RealVectorValue>("velocity")),
    
    // Initial values of the coordinates of the material deposition source
    _init_x_coords(getParam<std::vector<Real>>("init_x_coords")),
    _init_y_coords(getParam<std::vector<Real>>("init_y_coords")),
    _init_z_coords(getParam<std::vector<Real>>("init_z_coords")),
    
    // Total length during one scan
    _scan_length(getParam<std::vector<Real>>("scan_length")),
    
    // Threshold value of the ellipsoid function that activates the level set
	_level_set_activation_threshold(getParam<Real>("level_set_activation_threshold")),
	
	// _t_scan tracks the simulation time at which a new
    // scan begins after the condition based on the postprocessor
    // changes the coordinates of the heat source
    // imported from a VelocityEllipsoidHeatSource material object
    _t_scan(getMaterialProperty<Real>("t_scan")),
    
    // Material deposition track index
    // imported from a VelocityEllipsoidHeatSource material object
    _n_track(getMaterialProperty<int>("n_track"))
{
}

Real
VelocityEllipsoidAux::computeValue()
{
  // value of the level set variable at the previous time step
  Real old_level_set = _u[_qp];
  
  // Set initial coordinates for this track
  Real x_coord = _init_x_coords[_n_track[_qp]];
  Real y_coord = _init_y_coords[_n_track[_qp]];
  Real z_coord = _init_z_coords[_n_track[_qp]];
  
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);
  
  // center of the heat source
  Real x_t = x_coord + _velocity(0) * (_t - _t_scan[_qp]);
  Real y_t = y_coord + _velocity(1) * (_t - _t_scan[_qp]);
  Real z_t = z_coord + _velocity(2) * (_t - _t_scan[_qp]);
  
  // Calculate distance travelled by the heat source during this scan
  Real distance = std::pow(x_t - x_coord, 2);
  distance += std::pow(y_t - y_coord, 2);
  distance += std::pow(z_t - z_coord, 2);
  distance = std::sqrt(distance);
  
  // ellipsoid function value
  Real val;
  
  if (distance > _scan_length[_n_track[_qp]]) { // This single scan is over
	  
    val = 0.0;	  
	  
  } else {
	  
    val = 6.0 * std::sqrt(3.0) /
          (_rx * _ry * _rz * std::pow(libMesh::pi, 1.5)) *
          std::exp(-(3.0 * std::pow(x - x_t, 2.0) / std::pow(_rx, 2.0) +
                     3.0 * std::pow(y - y_t, 2.0) / std::pow(_ry, 2.0) +
                     3.0 * std::pow(z - z_t, 2.0) / std::pow(_rz, 2.0)));
  }
  
  if (val > _level_set_activation_threshold) { // ellipsoid function activating this _qp
	  
	  return _high_level_set_var;
	  
  } else {
	  
    if (old_level_set > _low_level_set_var) { // this was already activated
		 
      return _high_level_set_var;
		 
    } else { // this remains inactive

      return _low_level_set_var;

    }	 
  }
}
