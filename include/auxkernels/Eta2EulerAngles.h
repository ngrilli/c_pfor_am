// Nicolò Grilli
// Università di Bristol
// Fernando Valiente Dies
// ANSTO
// 26 Ottobre 2024

#pragma once

#include "AuxKernel.h"
#include "DelimitedFileReader.h"

/**
 * Output euler angles from eta_i phase fields
 */
class Eta2EulerAngles : public AuxKernel
{
public:
  static InputParameters validParams();

  Eta2EulerAngles(const InputParameters & parameters);

protected:
  virtual void assignEulerAngles();
  virtual Real computeValue();
  virtual void precalculateValue();

  /// Object providing the Euler angles
  ///const EulerAngleProvider & _euler;

  /// Grain tracker object
  ///const GrainTracker & _grain_tracker;

  /// Choice between phi1 Phi phi2
  MooseEnum _output_euler_angle;
  
  /// Euler angles file (in degrees) 
  /// each row must contain three Euler angles 
  /// which correspond to each grain orientation
  std::string _Euler_angles_file_name;
  
  /// Number of phase fields
  const unsigned int _op_num;

  /// Values of the phase field variables
  const std::vector<const VariableValue *> _vals;

  /// Euler angles for each phase field
  std::vector<RealVectorValue> _Euler_angles;

  /// precalculated element value
  Real _value;
};
