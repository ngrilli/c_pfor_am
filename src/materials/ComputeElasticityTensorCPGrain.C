// Nicolo Grilli
// ADSC Singapore
// 3 October 2020

#include "ComputeElasticityTensorCPGrain.h"
#include "RotationTensor.h"

registerMooseObject("TensorMechanicsApp", ComputeElasticityTensorCPGrain);

InputParameters
ComputeElasticityTensorCPGrain::validParams()
{
  InputParameters params = ComputeElasticityTensor::validParams();
  params.addClassDescription("Compute an elasticity tensor for crystal plasticity. "
                             "Euler angles are read from Euler angles input file "
							 "and can be assigned to physical volumes in GMSH");
  params.addParam<UserObjectName>("read_prop_user_object",
                                  "The GrainPropertyReadFile "
                                  "GeneralUserObject to read element "
                                  "specific property values from file");
  params.addCoupledVar("temp", 303.0,"Temperature, initialize at room temperature");
  params.addParam<Real>("reference_temperature",303.0,"reference temperature for elastic constants");
  params.addParam<Real>("dC11_dT", 0.0, "Change of C11 stiffness tensor component with temperature. "
									    "Use positive values, minus sign is added in the code. ");  
  params.addParam<Real>("dC12_dT", 0.0, "Change of C12 stiffness tensor component with temperature. "
								        "Use positive values, minus sign is added in the code. ");  
  params.addParam<Real>("dC44_dT", 0.0, "Change of C44 stiffness tensor component with temperature. "
									    "Use positive values, minus sign is added in the code. ");  
  return params;
}

ComputeElasticityTensorCPGrain::ComputeElasticityTensorCPGrain(const InputParameters & parameters)
  : ComputeElasticityTensor(parameters),
    _read_prop_user_object(isParamValid("read_prop_user_object")
                               ? &getUserObject<GrainPropertyReadFile>("read_prop_user_object")
                               : nullptr),
    _Euler_angles_mat_prop(declareProperty<RealVectorValue>("Euler_angles")),
	_temp(coupledValue("temp")),
    _reference_temperature(getParam<Real>("reference_temperature")),
	_dC11_dT(getParam<Real>("dC11_dT")),
	_dC12_dT(getParam<Real>("dC12_dT")),
	_dC44_dT(getParam<Real>("dC44_dT")),
    _crysrot(declareProperty<RankTwoTensor>("crysrot")),
    _R(_Euler_angles)
{
  // the base class guarantees constant in time, but in this derived class the
  // tensor will rotate over time once plastic deformation sets in
  revokeGuarantee(_elasticity_tensor_name, Guarantee::CONSTANT_IN_TIME);

  // the base class performs a passive rotation, but the crystal plasticity
  // materials use active rotation: recover unrotated _Cijkl here
  _Cijkl.rotate(_R.transpose());
}

void
ComputeElasticityTensorCPGrain::assignEulerAngles()
{
	
  if (_read_prop_user_object)
  {
    _Euler_angles_mat_prop[_qp](0) = _read_prop_user_object->getData(_current_elem, 0);
    _Euler_angles_mat_prop[_qp](1) = _read_prop_user_object->getData(_current_elem, 1);
    _Euler_angles_mat_prop[_qp](2) = _read_prop_user_object->getData(_current_elem, 2);
  }
  else
    _Euler_angles_mat_prop[_qp] = _Euler_angles;
  
}

void
ComputeElasticityTensorCPGrain::computeQpElasticityTensor()
{
  Real temp = _temp[_qp]; // Temperature
  Real deltatemp;
  
  // Properties assigned at the beginning of every call to material calculation
  assignEulerAngles();

  _R.update(_Euler_angles_mat_prop[_qp]);

  _crysrot[_qp] = _R.transpose();
  
  // Apply temperature dependence on _Cijkl
  // and save results on _Temp_Cijkl
  deltatemp = temp - _reference_temperature;
  temperatureDependence(deltatemp);
  
  _elasticity_tensor[_qp] = _Temp_Cijkl;

  _elasticity_tensor[_qp].rotate(_crysrot[_qp]);
}

// Temperature dependence of the elasticity tensor
// always referred to room temperature
void
ComputeElasticityTensorCPGrain::temperatureDependence(Real deltatemp)
{
  // Components with C11 coefficient
  _Temp_Cijkl(0, 0, 0, 0) = (1.0 - _dC11_dT * deltatemp) * _Cijkl(0, 0, 0, 0); // C1111
  _Temp_Cijkl(1, 1, 1, 1) = (1.0 - _dC11_dT * deltatemp) * _Cijkl(1, 1, 1, 1); // C2222
  _Temp_Cijkl(2, 2, 2, 2) = (1.0 - _dC11_dT * deltatemp) * _Cijkl(2, 2, 2, 2); // C3333

  // Components with C12 coefficient
  _Temp_Cijkl(0, 0, 1, 1) = (1.0 - _dC12_dT * deltatemp) * _Cijkl(0, 0, 1, 1); // C1122
  _Temp_Cijkl(1, 1, 0, 0) = (1.0 - _dC12_dT * deltatemp) * _Cijkl(1, 1, 0, 0);

  _Temp_Cijkl(0, 0, 2, 2) = (1.0 - _dC12_dT * deltatemp) * _Cijkl(0, 0, 2, 2); // C1133
  _Temp_Cijkl(2, 2, 0, 0) = (1.0 - _dC12_dT * deltatemp) * _Cijkl(2, 2, 0, 0);

  _Temp_Cijkl(1, 1, 2, 2) = (1.0 - _dC12_dT * deltatemp) * _Cijkl(1, 1, 2, 2); // C2233
  _Temp_Cijkl(2, 2, 1, 1) = (1.0 - _dC12_dT * deltatemp) * _Cijkl(2, 2, 1, 1);

  // Components with C44 coefficient
  _Temp_Cijkl(1, 2, 1, 2) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(1, 2, 1, 2); // C2323
  _Temp_Cijkl(2, 1, 2, 1) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(2, 1, 2, 1);
  _Temp_Cijkl(2, 1, 1, 2) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(2, 1, 1, 2);
  _Temp_Cijkl(1, 2, 2, 1) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(1, 2, 2, 1);

  _Temp_Cijkl(0, 2, 0, 2) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(0, 2, 0, 2); // C1313
  _Temp_Cijkl(2, 0, 2, 0) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(2, 0, 2, 0);
  _Temp_Cijkl(2, 0, 0, 2) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(2, 0, 0, 2);
  _Temp_Cijkl(0, 2, 2, 0) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(0, 2, 2, 0);

  _Temp_Cijkl(0, 1, 0, 1) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(0, 1, 0, 1); // C1212
  _Temp_Cijkl(1, 0, 1, 0) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(1, 0, 1, 0);
  _Temp_Cijkl(1, 0, 0, 1) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(1, 0, 0, 1);
  _Temp_Cijkl(0, 1, 1, 0) = (1.0 - _dC44_dT * deltatemp) * _Cijkl(0, 1, 1, 0); 
  
}

