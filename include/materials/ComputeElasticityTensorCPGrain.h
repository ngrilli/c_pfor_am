// Nicolo Grilli
// ADSC Singapore
// 3 October 2020

#pragma once

#include "ComputeElasticityTensor.h"
#include "GrainPropertyReadFile.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"

/**
 * ComputeElasticityTensorCPGrain defines an elasticity tensor material object for crystal plasticity.
 * It is based on the user object GrainPropertyReadFile
 * that allows to assign grains from the euler angles input file
 * to the physical volumes in GMSH
 */
class ComputeElasticityTensorCPGrain : public ComputeElasticityTensor
{
public:
  static InputParameters validParams();

  ComputeElasticityTensorCPGrain(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor() override;

  virtual void assignEulerAngles();
  
  virtual void temperatureDependence(Real);

  /**
   * Element property read user object
   * Presently used to read Euler angles -  see test
   */
  const GrainPropertyReadFile * const _read_prop_user_object;
  
  MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;
  
  const VariableValue & _temp;
  
  const Real _reference_temperature;
  
  const Real _dC11_dT;
  const Real _dC12_dT;
  const Real _dC44_dT;

  /// Crystal Rotation Matrix
  MaterialProperty<RankTwoTensor> & _crysrot;

  /// Rotation matrix
  RotationTensor _R;
  
  /// Stiffness tensor modified by temperature
  RankFourTensor _Temp_Cijkl;
};
