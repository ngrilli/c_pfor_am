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
  Real TempValue = 0.0;
  Real x, y, z;
  Real x0, y0, z0;

  Real LaserWidth = 200.0; // characteristic diameter of laser (micrometres)
  Real WidthBehind = 400.0; // characteristic length scale behind of the laser tail (micrometres)
  Real MaxTemp = 779.0; // max temperature above 293.0 in the tail of the laser
  Real LaserShrinkBegin = 390.0; // beginning of laser shrinking with respect to laser tail (micrometres)
  Real TailBegin = 740.0; // beginning of the laser tail (micrometres)
  Real LaserCentre = 790.0; // middle of the laser tail, where laser starts to shrink (micrometres)
  Real LaserMaxT = 940.0; // End of the laser tail
  Real LaserHead = 1040.0;
  Real zLength = 90.0; // characteristic length of high temperature region in depth (micrometres)
  Real HeatingTime = 1.0e-4; // Characteristic time for heating up (s)
  
  // Get coordinates of the current IP
  x = _q_point[_qp](0);
  y = _q_point[_qp](1);
  z = _q_point[_qp](2);
  
  if (_laser_init_coord.size() != 3)
    mooseError("LaserScanAux: Error in reading laser_init_coord: "
               "must be in the format 'x0 y0 z0', "
			   "where z0 is on the upper surface. ");
  
  // Calculate laser centre as a function of time
  x0 = _laser_init_coord[0] + _scan_velocity * _t;
  y0 = _laser_init_coord[1];
  z0 = _laser_init_coord[2];
  LaserShrinkBegin += x0;
  TailBegin += x0;
  LaserCentre += x0;
  LaserMaxT += x0;
  LaserHead += x0;
  
  // different sections of the temperature profile
  // depending on coordinates
  if (x < x0) {
    LaserWidth *= 2.0;
	TempValue += MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0)) * 
	                       std::exp((x-x0)/WidthBehind);
} else if (x >= x0 && x < LaserShrinkBegin) {
    LaserWidth *= 2.0;
	TempValue += MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0));
} else if (x >= LaserShrinkBegin && x < TailBegin) {	
    LaserWidth *= (2.0 - 1.5 * std::pow((x-LaserShrinkBegin)/(TailBegin-LaserShrinkBegin),0.3));
    TempValue += MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0));
} else if (x >= TailBegin && x < LaserCentre) {
	LaserWidth *= 0.4;
	TempValue += MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0));
} else if (x >= LaserCentre && x < LaserMaxT) {
	LaserWidth *= 0.4 * (1.0 - 0.4 * (x-LaserCentre)/(LaserMaxT-LaserCentre));
	MaxTemp *= (1.0 + 0.25 * (x-LaserCentre)/(LaserMaxT-LaserCentre));
	TempValue += MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0));
} else if (x >= LaserMaxT && x < LaserHead) {
	LaserWidth *= 0.4 * 0.6 * (1.0 - std::pow(std::max((x-LaserMaxT)/(LaserHead-LaserMaxT),1.0e-40),0.2));
    MaxTemp *= 1.25 * (1.0 - 0.5 * (x-LaserMaxT)/(LaserHead-LaserMaxT));
	TempValue += MaxTemp * std::exp(-std::pow((y-y0)/LaserWidth,2.0));
} else {
	TempValue += 0.0;
}

  // Exponential decrease of the temperature in depth
  TempValue *= std::exp(-(std::abs(z-z0))/zLength);
  
  if (_t < HeatingTime) {
	TempValue *= (_t / HeatingTime); 
  }
  
  // Add room temperature
  TempValue += 293.0;

  return TempValue;
}
