// Nicolò Grilli
// Università di Bristol
// 23 Novembre 2024

#pragma once

#include "MeshGenerator.h"
#include "DelimitedFileReader.h"

/**
 * MeshGenerator for assigning subdomain IDs of all elements
 * by reading the list of IDs for all elements from file
 */
class SubdomainFileReader : public MeshGenerator
{
public:
  static InputParameters validParams();

  SubdomainFileReader(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// mesh to modify
  std::unique_ptr<MeshBase> & _input;
  
  /// File should contain the subdomain ID for each element
  std::string _subdomain_id_file_name;
};
