// Nicolò Grilli
// 1 Febbraio 2021
// National University of Singapore

#include "DislocationLoopsIC.h"
#include "libmesh/utility.h"

#include <cmath>
#include <fstream>

registerMooseObject("MooseApp", DislocationLoopsIC);

defineLegacyParams(DislocationLoopsIC);

InputParameters
DislocationLoopsIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addParam<std::vector<Real>>("centrex", "X coordinates of the centres of the loops");
  params.addParam<std::vector<Real>>("centrey", "Y coordinates of the centres of the loops");
  params.addParam<std::vector<Real>>("radii", "Radii of the loops, signed for positive/negative loops");
  params.addParam<std::vector<Real>>("width", "Width of the loops");
  params.addParam<std::vector<Real>>("rho_max", "Max dislocation density in one loop");
  params.addParam<std::vector<Real>>("depth", "Depth at which slip planes are placed");
  params.addParam<std::vector<Real>>("thickness", "Thickness of the slip planes");
  MooseEnum var_options("rhotot rhoedgegnd rhoscrewgnd qtot", "rhotot");
  params.addParam<MooseEnum>("variable_type",
                             var_options,
                             "Type of variable for the dislocation loop."); 
  params.addParam<int>("slip_sys_index", 0, "Slip system index to determine slip direction "
							   "for instance from 0 to 11 for FCC.");
  params.addParam<UserObjectName>("read_prop_user_object",
                                  "The GrainPropertyReadFile "
                                  "GeneralUserObject to read element "
                                  "specific property values from file");
  params.addRequiredParam<FileName>("slip_sys_file_name",
                                    "Name of the file containing the slip system");		
  params.addRequiredParam<int>("nss", "Number of slip systems");	
  params.addParam<bool>("is3D",false,"Initial dislocation density on parallel slip planes in 3D");  
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
	_depth(getParam<std::vector<Real>>("depth")),
	_thickness(getParam<std::vector<Real>>("thickness")),
    _variable_type(getParam<MooseEnum>("variable_type")),
	_slip_sys_index(getParam<int>("slip_sys_index")),
    _read_prop_user_object(isParamValid("read_prop_user_object")
                               ? &getUserObject<GrainPropertyReadFile>("read_prop_user_object")
                               : nullptr),
    _slip_sys_file_name(getParam<FileName>("slip_sys_file_name")),
    _nss(getParam<int>("nss")),
	_is3D(getParam<bool>("is3D")),
    _mo(_nss * LIBMESH_DIM),
    _no(_nss * LIBMESH_DIM),
	_rot_mo(_nss * LIBMESH_DIM),
	_rot_no(_nss * LIBMESH_DIM),
	_rot_to(_nss * LIBMESH_DIM),
	_R(_Euler_angles_sc)
{
}

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
  Real absradius; // Absolute value of the radius of curvature 
  Real val;
  
  bool OutOfSlipPlane = true; // point is out of the slip plane -> no initialisation
  
  Point pp;	
	
  if (_centrex.size() <= 0 || _centrey.size() <= 0 || _radii.size() <= 0 || _width.size() <= 0)
    mooseError("Error in reading dislocation loops properties");

  if (_is3D) {
    if (_depth.size() <= 0 || _thickness.size() <= 0) {
	  mooseError("Error in reading parallel slip planes properties");	
	}	  
	if (_depth.size() != _thickness.size()) {
	  mooseError("Size of _depth must be the same as size of _thickness");	
	}
  }

  // get coordinates on the slip system
  pp = projectOnSlipPlane(p);

  for (unsigned int i = 0; i < _centrex.size(); ++i) {
	  
	absradius = std::abs(_radii[i]);
	  
	R = std::sqrt((pp(0)-_centrex[i])*(pp(0)-_centrex[i])+(pp(1)-_centrey[i])*(pp(1)-_centrey[i]));
	
	rhototloop = std::exp(-(R-absradius)*(R-absradius)/(2.0*_width[i]*_width[i]));
	rhototloop *= (1.0/(2.506627 * _width[i]));
	rhototloop *= _rho_max[i];
	
	rhotot += rhototloop; 
	qtot = qtot + (rhototloop / _radii[i]);
	
	// Multiply total dislocation density by the unit radial vector
	// to obtain dislocation density vector
	if (R > 0.0001) {
	  radialvx = ((pp(0)-_centrex[i]) / R) * std::copysign(1.0, _radii[i]);
	  radialvy = ((pp(1)-_centrey[i]) / R) * std::copysign(1.0, _radii[i]);		
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
  
  // catch possible NaN values
  if (val != val) {
	mooseWarning("Projected coordinate X ", pp(0));
	mooseWarning("Projected coordinate Y ", pp(1));
	mooseWarning("Projected coordinate Z ", pp(2));
	mooseWarning("Distance from point R ", R);
	mooseWarning("Total dislocation density ", rhotot);
	mooseWarning("Edge GND density ", rhoedgegnd);
	mooseWarning("Screw GND density ", rhoscrewgnd);
	mooseWarning("Curvature density ", qtot);
	mooseError("NaN detected in DislocationLoopsIC");
  }
  
  // return 0 if the point is out of the selected
  // parallel slip planes in 3D
  if (_is3D) {
    for (unsigned int i = 0; i < _depth.size(); ++i) {
      if ((pp(2) > _depth[i]) && (pp(2) < (_depth[i] + _thickness[i]))) {
	    OutOfSlipPlane = false;	  
      }	
	}	  
  } else {
	OutOfSlipPlane = false;  
  }
  
  if (OutOfSlipPlane) {
	  
    return 0.0;
	
  } else {

    return val; 

  }

}

// project an arbitrary point p on the slip plane
// passing through that point
Point
DislocationLoopsIC::projectOnSlipPlane(const Point & p)
{
  Point projectedp;	

  assignEulerAngles();
  
  getSlipSystems();
  
  rotateSlipSystems();
	
  projectedp(0) = p(0) * _rot_mo(_slip_sys_index * LIBMESH_DIM);
  projectedp(0) += p(1) * _rot_mo(_slip_sys_index * LIBMESH_DIM + 1);
  projectedp(0) += p(2) * _rot_mo(_slip_sys_index * LIBMESH_DIM + 2);
  
  // _rot_to points towards -y direction
  // therefore a change of sign is needed
  projectedp(1) = p(0) * (-1.0) * _rot_to(_slip_sys_index * LIBMESH_DIM);
  projectedp(1) += p(1) * (-1.0) * _rot_to(_slip_sys_index * LIBMESH_DIM + 1);
  projectedp(1) += p(2) * (-1.0) * _rot_to(_slip_sys_index * LIBMESH_DIM + 2);
  
  projectedp(2) = p(0) * _rot_no(_slip_sys_index * LIBMESH_DIM); 
  projectedp(2) += p(1) * _rot_no(_slip_sys_index * LIBMESH_DIM + 1);
  projectedp(2) += p(2) * _rot_no(_slip_sys_index * LIBMESH_DIM + 2);
	
  return projectedp;
}

void
DislocationLoopsIC::assignEulerAngles()
{
  if (_read_prop_user_object)
  {
    _Euler_angles_sc(0) = _read_prop_user_object->getData(_current_elem, 0);
    _Euler_angles_sc(1) = _read_prop_user_object->getData(_current_elem, 1);
    _Euler_angles_sc(2) = _read_prop_user_object->getData(_current_elem, 2);
  }
  else
    mooseError("DislocationLoopsIC: Error in reading Euler angles");
}

// Read slip systems from file
void
DislocationLoopsIC::getSlipSystems()
{
  Real vec[LIBMESH_DIM];
  std::ifstream fileslipsys;
  
  MooseUtils::checkFileReadable(_slip_sys_file_name);

  fileslipsys.open(_slip_sys_file_name.c_str());

  for (unsigned int i = 0; i < _nss; ++i)
  {
    // Read the slip normal
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      if (!(fileslipsys >> vec[j]))
        mooseError("Crystal Plasticity Error: Premature end of file reading slip system file \n");

    // Normalize the vectors
    Real mag;
    mag = Utility::pow<2>(vec[0]) + Utility::pow<2>(vec[1]) + Utility::pow<2>(vec[2]);
    mag = std::sqrt(mag);

    for (unsigned j = 0; j < LIBMESH_DIM; ++j)
      _no(i * LIBMESH_DIM + j) = vec[j] / mag;

    // Read the slip direction
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      if (!(fileslipsys >> vec[j]))
        mooseError("Crystal Plasticity Error: Premature end of file reading slip system file \n");

    // Normalize the vectors
    mag = Utility::pow<2>(vec[0]) + Utility::pow<2>(vec[1]) + Utility::pow<2>(vec[2]);
    mag = std::sqrt(mag);

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      _mo(i * LIBMESH_DIM + j) = vec[j] / mag;

    mag = 0.0;
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      mag += _mo(i * LIBMESH_DIM + j) * _no(i * LIBMESH_DIM + j);

    if (std::abs(mag) > 1e-8)
      mooseError(
          "Crystal Plasicity Error: Slip direction and normal not orthonormal, System number = ",
          i,
          "\n");
  }

  fileslipsys.close();  
}

// rotate slip systems directions and normals 
// from the lattice frame
// to the sample frame
void
DislocationLoopsIC::rotateSlipSystems()
{
  RealVectorValue temp_mo;
  RealVectorValue temp_no;
  RealVectorValue temp_to;
  
  // update function builds the passive rotation
  _R.update(_Euler_angles_sc);  
  _Rt = _R.transpose();

  for (unsigned int i = 0; i < _nss; ++i) {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j) {
	  _rot_mo(i * LIBMESH_DIM + j) = 0.0;
	  _rot_no(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k) {
        _rot_mo(i * LIBMESH_DIM + j) += _Rt(j, k) * _mo(i * LIBMESH_DIM + k);
		_rot_no(i * LIBMESH_DIM + j) += _Rt(j, k) * _no(i * LIBMESH_DIM + k);
      }	  
    }
  }
  
  // Store slip direction for screw dislocations
  // already normalized
  for (unsigned int i = 0; i < _nss; ++i)
  {
	for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
	{
	  temp_mo(j) = _rot_mo(i * LIBMESH_DIM + j);
	  temp_no(j) = _rot_no(i * LIBMESH_DIM + j);
	}		
	
	temp_to = temp_mo.cross(temp_no);
	  
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
  	  _rot_to(i * LIBMESH_DIM + j) = temp_to(j);
  	}
  }
  
}

