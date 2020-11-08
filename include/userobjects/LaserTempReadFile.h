// Nicolo Grilli
// National University of Singapore
// 8 Novembre 2020

#pragma once

#include "GeneralUserObject.h"

/**
 * Read temperature field from file
 * Usable for generated mesh
 */

class LaserTempReadFile : public GeneralUserObject
{
public:
  static InputParameters validParams();

  LaserTempReadFile(const InputParameters & parameters);
  virtual ~LaserTempReadFile() {}

  virtual void initialize() {}
  virtual void execute() {}
  virtual void finalize() {}

  /**
   * This function reads element data from file
   */
  void readElementData();

  /**
   * This function assign temperature data to elements
   */
  Real getData(const Elem *, unsigned int) const;

protected:
  
  /// Name of file containing temperature values
  std::string _temperature_file_name;
  
  /// Number of temperature data field in time
  const unsigned int _temperature_num_step;
  
  MooseMesh & _mesh;
  
  /// Store temperature values read from file
  std::vector<Real> _data;

private:
  unsigned int _nelem;
  Point _top_right;
  Point _bottom_left;
  Point _range;
  Real _max_range;
};
