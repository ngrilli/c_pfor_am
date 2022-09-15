// Nicolò Grilli
// University of Bristol
// 4 Settembre 2022

#pragma once

#include "Material.h"
#include "ElementPropertyReadFile.h"

/**
 * Parse a scalar material property from file
 * and assign to elements based on element number
 */
class FileParsedMaterial : public Material
{
public:
  static InputParameters validParams();

  FileParsedMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  
  /// arbitrary name of the material property
  const std::string _prop_name;

  /// Misorientation between neighboring grains
  MaterialProperty<Real> & _property;
  
  /// property read user object used to read the material property
  /// File must contain one scalar in each row
  /// corresponding to the property in the element with index
  /// equal to the row number
  const ElementPropertyReadFile * const _read_prop_user_object;
};
