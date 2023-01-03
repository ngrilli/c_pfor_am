// Fernando Valiente Dies & Nicolò Grilli
// University of Sydney & University of Bristol
// 26 Novembre 2022

#include "VelocityEllipsoidHeatSource.h"

registerMooseObject("HeatConductionApp", VelocityEllipsoidHeatSource);

InputParameters
VelocityEllipsoidHeatSource::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Double ellipsoid volumetric source heat, the motion is determined "
                             "by input velocity, starting positions and a postprocessor. ");
  params.addRequiredParam<Real>("power", "power");
  params.addParam<Real>("efficiency", 1, "process efficiency");
  params.addRequiredParam<Real>("rx", "effective transverse ellipsoid radius");
  params.addRequiredParam<Real>("ry", "effective longitudinal ellipsoid radius");
  params.addRequiredParam<Real>("rz", "effective depth ellipsoid radius");
  params.addParam<Real>(
      "factor", 1, "scaling factor that is multiplied to the heat source to adjust the intensity");
  params.addRequiredParam<RealVectorValue>("velocity", "Velocity vector");
  
  // Every time the postprocessor condition is satisfied, the heat source is moved to the next set of coordinates
  params.addRequiredParam<std::vector<Real>>("init_x_coords", "Initial values of x coordinates of the heat source");
  params.addRequiredParam<std::vector<Real>>("init_y_coords", "Initial values of y coordinates of the heat source");
  params.addRequiredParam<std::vector<Real>>("init_z_coords", "Initial values of z coordinates of the heat source");
  
  params.addRequiredParam<PostprocessorName>("temperature_pp","Postprocessor with temperature value to determine heat source motion.");
  
  params.addRequiredParam<std::vector<Real>>("scan_length","Total length during one scan. "
                                                           "After this length the laser is switched off. ");

  params.addRequiredParam<Real>("threshold_temperature","When the temperature provided by the postprocessor decreases "
                                                        "below this threshold, the heat source is moved to the next "
                                                        "set of coordinates. ");                                                   
  return params;
}

VelocityEllipsoidHeatSource::VelocityEllipsoidHeatSource(const InputParameters & parameters)
  : Material(parameters),
    _P(getParam<Real>("power")),
    _eta(getParam<Real>("efficiency")),
    _rx(getParam<Real>("rx")),
    _ry(getParam<Real>("ry")),
    _rz(getParam<Real>("rz")),
    _f(getParam<Real>("factor")),
    _velocity(getParam<RealVectorValue>("velocity")), // Scanning speed vector
    
    // Initial values of the coordinates of the heat source
    _init_x_coords(getParam<std::vector<Real>>("init_x_coords")),
    _init_y_coords(getParam<std::vector<Real>>("init_y_coords")),
    _init_z_coords(getParam<std::vector<Real>>("init_z_coords")),
    
    // Postprocess with temperature value
    _temperature_pp(getPostprocessorValue("temperature_pp")),

    // _t_scan tracks the simulation time at which a new
    // scan begins after the condition based on the postprocessor
    // changes the coordinates of the heat source
    _t_scan(declareProperty<Real>("t_scan")),
    
    // Total time during one scan
    _scan_length(getParam<std::vector<Real>>("scan_length")),
    
    // Threshold temperature for the postprocessor condition
    _threshold_temperature(getParam<Real>("threshold_temperature")),
    
    // Heat source track index
    _n_track(declareProperty<int>("n_track")),
    
    // Volumetric heat source used by the kernel
    _volumetric_heat(declareADProperty<Real>("volumetric_heat"))
{
}

void
VelocityEllipsoidHeatSource::initQpStatefulProperties()
{
  // Initialize time tracking and number of tracks
  _t_scan[_qp] = _t;
  _n_track[_qp] = 0;
}

void
VelocityEllipsoidHeatSource::computeQpProperties()
{	
  // When maximum number of tracks is reached
  // heat source is switched off
  if (_n_track[_qp] >= _init_x_coords.size()) {
	  
    _volumetric_heat[_qp] = 0.0;	
	return;
	
  }

  // Set initial coordinates for this track
  _x_coord = _init_x_coords[_n_track[_qp]];
  _y_coord = _init_y_coords[_n_track[_qp]];
  _z_coord = _init_z_coords[_n_track[_qp]];

  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);
  
  // center of the heat source
  Real x_t = _x_coord + _velocity(0) * (_t - _t_scan[_qp]);
  Real y_t = _y_coord + _velocity(1) * (_t - _t_scan[_qp]);
  Real z_t = _z_coord + _velocity(2) * (_t - _t_scan[_qp]);

  // Calculate distance travelled by the heat source during this scan
  Real distance = std::pow(x_t - _x_coord, 2);
  distance += std::pow(y_t - _y_coord, 2);
  distance += std::pow(z_t - _z_coord, 2);
  distance = std::sqrt(distance);

  if (distance > _scan_length[_n_track[_qp]]) { // This single scan is over
	  
    _volumetric_heat[_qp] = 0.0;	
    checkPPcondition();  
	  
  } else {

    _volumetric_heat[_qp] = 6.0 * std::sqrt(3.0) * _P * _eta * _f /
                            (_rx * _ry * _rz * std::pow(libMesh::pi, 1.5)) *
                            std::exp(-(3.0 * std::pow(x - x_t, 2.0) / std::pow(_rx, 2.0) +
                                       3.0 * std::pow(y - y_t, 2.0) / std::pow(_ry, 2.0) +
                                       3.0 * std::pow(z - z_t, 2.0) / std::pow(_rz, 2.0)));
  }

}

// Check if the postprocessor temperature condition is satisfied
// and change the initial coordinates and scan time
void
VelocityEllipsoidHeatSource::checkPPcondition()
{
  // cooling condition currently not implemented
  // There is no condition to check that the temperature is decreasing
  // when the laser position is changed, only threshold is used
  
  if (_temperature_pp < _threshold_temperature) { // reached threshold temperature
		
    // update initial heat source coordinate and track time
    _n_track[_qp] += 1;
    _t_scan[_qp] = _t;

  }
}
