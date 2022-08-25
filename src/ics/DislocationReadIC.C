// Daijun Hu
// National University of Singapore
// Nicolò Grilli
// Univerisity of Bristol
// 25 Agosto 2022

#include "DislocationReadIC.h"
#include "libmesh/utility.h"

#include <cmath>
#include <fstream>

registerMooseObject("MooseApp", DislocationReadIC);

InputParameters
DislocationReadIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addParam<UserObjectName>("read_dislocation_user_object",
                                  "The PropertyReadFile "
                                  "GeneralUserObject to read element "
                                  "specific dislocation values from file. "
								  "File must contain four columns with: "
								  "rho_tot rho_e rho_s q ");  
  MooseEnum var_options("rhotot rhoedgegnd rhoscrewgnd qtot", "rhotot");
  params.addParam<MooseEnum>("variable_type",
                             var_options,
                             "Type of variable for the dislocation initialization."); 
  params.addClassDescription("Initialize dislocation densities and curvature read from file. ");
  return params;
}

DislocationReadIC::DislocationReadIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _read_dislocation_user_object(isParamValid("read_dislocation_user_object")
                               ? &getUserObject<PropertyReadFile>("read_dislocation_user_object")
                               : nullptr),
    _variable_type(getParam<MooseEnum>("variable_type"))
{
}

Real
DislocationReadIC::value(const Point & p)
{
  Real val; // output variable
  
  // temporary variables to store the values read from file
  Real rhotot;
  Real rhoedgegnd;
  Real rhoscrewgnd;
  Real qtot;
  
  if (_read_dislocation_user_object)
  {
    rhotot = _read_dislocation_user_object->getData(_current_elem, 0);
    rhoedgegnd = _read_dislocation_user_object->getData(_current_elem, 1);
    rhoscrewgnd = _read_dislocation_user_object->getData(_current_elem, 2);
    qtot = _read_dislocation_user_object->getData(_current_elem, 3);
  }
  else
    mooseError("DislocationLoopsIC: Error in reading Euler angles");

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
	mooseError("NaN detected in DislocationReadIC");
  }

  return val;
}
