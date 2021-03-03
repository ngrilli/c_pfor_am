// Nicolò Grilli
// 1 Febbraio 2021
// National University of Singapore

#pragma once

#include "InitialCondition.h"
#include "GrainPropertyReadFile.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"

// Forward Declarations
class DislocationLoopsIC;
class InputParameters;

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<DislocationLoopsIC>();

/**
 * Initialize dislocation loops.
 */
class DislocationLoopsIC : public InitialCondition
{
public:
  static InputParameters validParams();

  DislocationLoopsIC(const InputParameters & parameters);

protected:

  virtual Real value(const Point & p) override;
  
  virtual Point projectOnSlipPlane(const Point & p);
  
  virtual void assignEulerAngles();
  
  virtual void getSlipSystems();
  
  virtual void rotateSlipSystems();
  
  std::vector<Real> _centrex;
  std::vector<Real> _centrey;
  std::vector<Real> _radii; // signed for positive/negative loops
  std::vector<Real> _width;
  std::vector<Real> _rho_max;
  std::vector<Real> _depth;
  std::vector<Real> _thickness;
  
  // Type of variable for the dislocation loop
  MooseEnum _variable_type;
  
  // Slip system index to determine slip direction
  const unsigned int _slip_sys_index;
  
  /**
   * Element property read user object
   * Presently used to read Euler angles
   */
  const GrainPropertyReadFile * const _read_prop_user_object;
  
  ///File should contain slip plane normal and direction
  /// not yet normalized
  std::string _slip_sys_file_name;
  
  /// Number of slip system
  const unsigned int _nss;
  
  /// Dislocation loops on parallel slip planes in 3D
  bool _is3D;  
  
  /// Euler angles
  /// works only for single crystal
  RealVectorValue _Euler_angles_sc;
  
  /// Normalized slip direction and normal and edge dislocation line
  /// in the crystal lattice frame
  DenseVector<Real> _mo;
  DenseVector<Real> _no;

  /// Normalized slip direction and normal and edge dislocation line
  /// in the sample frame  
  DenseVector<Real> _rot_mo;
  DenseVector<Real> _rot_no;
  DenseVector<Real> _rot_to;

  /// Rotation matrix (passive)
  RotationTensor _R;
  
  /// Rotation matrix (active) 
  RankTwoTensor _Rt;

};
