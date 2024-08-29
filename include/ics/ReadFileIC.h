// Nicolò Grilli
// Università di Bristol
// 29 Agosto 2024

#pragma once

#include "InitialCondition.h"
#include "DelimitedFileReader.h"

class InputParameters;

template <typename T>
InputParameters validParams();

/**
 * Defines a boundary condition that is read from file
 * and assigned to the nodes of a structured mesh
 */
class ReadFileIC : public InitialCondition
{
public:
  static InputParameters validParams();

  ReadFileIC(const InputParameters & parameters);

protected:

  /**
   * The value of the variable at a point.
   */
  virtual Real value(const Point & p) override;

  /**
   * The value of the gradient at a point.
   */
  virtual RealGradient gradient(const Point & p) override;

  /// File should contain a single column with the value
  /// of the variable in each element of a structured mesh
  std::string _ic_file_name;
  
  /// Element size in the structured mesh
  /// assuming size along x and y are the same
  const Real _element_size;
  
  /// Number of elements along the x and y axes in the structured mesh
  const unsigned int _nx;
  const unsigned int _ny;
};
