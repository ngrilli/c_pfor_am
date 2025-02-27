// Nicolò Grilli
// Università di Bristol
// Fernando Valiente Dies
// ANSTO
// 26 Febbraio 2025

#pragma once

#include "Material.h"
#include "DelimitedFileReader.h"

/**
 * Output euler angles from eta_i phase fields
 * into a material property
 */
class Eta2EulerAnglesMaterial : public Material
{
public:
  static InputParameters validParams();

  Eta2EulerAnglesMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  
  /// assign Euler angles read from file to material property
  virtual void assignEulerAngles();

  /// Euler angles file (in degrees) 
  /// each row must contain three Euler angles 
  /// which correspond to each grain orientation
  const std::string _Euler_angles_file_name;
  
  /// Number of phase fields
  const unsigned int _op_num;

  /// Phase field variables
  const std::vector<const VariableValue *> _vals;
  
  /// Euler angles for each phase field
  std::vector<RealVectorValue> _Euler_angles;
  
  /// Material property that stores the values of the Euler Angles for postprocessing
  MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;
};
