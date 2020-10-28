// Nicolò Grilli
// 28 Ottobre 2020
// National University of Singapore

#include "LaserScanAux.h"

registerMooseObject("TensorMechanicsApp", LaserScanAux);

InputParameters
LaserScanAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Compute temperature due to laser scan during SLM from CFD simulation");
  params.addParam<std::string>("base_name", "Mechanical property base name");
  params.addParam<Real>("scan_velocity",0.0,"Laser scan velocity magnitude, assumed along negative x");
  params.addRequiredParam<std::vector<Real>>("laser_init_coord","Initial (x0,y0,z0) coordinates of the laser. "
                                                                "z0 must be on the upper surface. ");
  return params;
}

LaserScanAux::LaserScanAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _scan_velocity(getParam<Real>("scan_velocity")),
    _laser_init_coord(getParam<std::vector<Real>>("laser_init_coord"))
{
}

// Calculate temperature based on fitting from CFD simulation
Real
LaserScanAux::computeValue()
{
  // Some of these parameters may eventually go in the input file
  Real TempValue = 293.0;
  Real x, y, z;
  Real x0, y0, z0;
  Real LaserWidth = 200.0; // characteristic diameter of laser at laser tail (micrometres)
  Real WidthAhead = 300.0; // characteristic length scale ahead of the laser centre (micrometres)
  Real MaxTemp = 755.0; // max temperature above 298.0 in the centre of the laser
  Real LaserShrinkBegin = 104.21; // beginning of laser shrinking with respect to laser centre (micrometres)
  Real TailBegin = 208.42; // beginning of the laser tail (micrometres)
  Real TailMiddle = 547.11; // middle of the laser tail, where laser starts to shrink (micrometres)
  Real TailEnd = 997.11; // End of the laser tail
  Real zLength = 90.0; // characteristic length of high temperature region in depth (micrometres)
  Real HeatingTime = 1.0; // Characteristic time for heating up (s)
  
  // Get coordinates of the current IP
  x = _q_point[_qp](0);
  y = _q_point[_qp](1);
  z = _q_point[_qp](2);
  
  if (_laser_init_coord.size() != 3)
    mooseError("LaserScanAux: Error in reading laser_init_coord: "
               "must be in the format 'x0 y0 z0', "
			   "where z0 is on the upper surface. ");
  
  // Calculate laser centre as a function of time
  x0 = _laser_init_coord[0] - _scan_velocity * _t;
  y0 = _laser_init_coord[1];
  z0 = _laser_init_coord[2];
  LaserShrinkBegin += x0;
  TailBegin += x0;
  TailMiddle += x0;
  TailEnd += x0;

  // different sections of the temperature profile
  // depending on coordinates
  if (x < x0) {
    LaserWidth *= 2.0;
	TempValue += MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0)) * 
	                       std::exp((x-x0)/WidthAhead);
} else if (x >= x0 && x < LaserShrinkBegin) {
    LaserWidth *= 2.0;
	TempValue += MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0));
} else if (x >= LaserShrinkBegin && x < TailBegin) {	
    LaserWidth *= (1.0 + (TailBegin-x)/(TailBegin-LaserShrinkBegin));
    TempValue += MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0));
} else if (x >= TailBegin && x < TailMiddle) {
	TempValue += 0.95 * MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0));
} else if (x >= TailMiddle && x < TailEnd) {
	LaserWidth *= std::max(1.0 - (x-TailMiddle)/(TailEnd-TailMiddle),0.001);
	TempValue += 0.95 * MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0)) * 
	             std::pow(1.0 - (x-TailMiddle)/(TailEnd-TailMiddle),0.2);
} else {
	TempValue += 0.0;
}

  // Exponential decrease of the temperature in depth
  TempValue *= std::exp(-(std::abs(z-z0))/zLength);
  
  if (_t < 1.0) {
	TempValue *= _t; 
  }
  
  return TempValue;
}
