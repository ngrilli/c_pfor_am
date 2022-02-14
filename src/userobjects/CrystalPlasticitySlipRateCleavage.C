// Nicolò Grilli
// University of Bristol
// 13 Febbraio 2022

// this is the same as CrystalPlasticitySlipRateGSS
// but can output the slip plane normal using the
// function calcFlowDirection

#include "CrystalPlasticitySlipRateCleavage.h"

#include <fstream>

registerMooseObject("TensorMechanicsApp", CrystalPlasticitySlipRateCleavage);

InputParameters
CrystalPlasticitySlipRateCleavage::validParams()
{
  InputParameters params = CrystalPlasticitySlipRateGSS::validParams();
  params.addClassDescription("Phenomenological constitutive model slip rate class. Override the "
                             "virtual functions in your class. It can output the slip plane " 
							 "normal using the function calcFlowDirection.");
  return params;
}

CrystalPlasticitySlipRateCleavage::CrystalPlasticitySlipRateCleavage(const InputParameters & parameters)
  : CrystalPlasticitySlipRateGSS(parameters)
{
  if (_slip_sys_flow_prop_file_name.length() != 0)
    readFileFlowRateParams();
  else
    getFlowRateParams();
}

void
CrystalPlasticitySlipRateCleavage::calcFlowDirection(unsigned int qp,
                                                std::vector<RankTwoTensor> & flow_direction,
												std::vector<Real> & slip_plane_normals) const
{
  DenseVector<Real> mo(LIBMESH_DIM * _variable_size), no(LIBMESH_DIM * _variable_size);

  // Update slip direction and normal with crystal orientation
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo(i * LIBMESH_DIM + j) =
            mo(i * LIBMESH_DIM + j) + _crysrot[qp](j, k) * _mo(i * LIBMESH_DIM + k);
    }

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        no(i * LIBMESH_DIM + j) =
            no(i * LIBMESH_DIM + j) + _crysrot[qp](j, k) * _no(i * LIBMESH_DIM + k);
    }
  }

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _variable_size; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        flow_direction[i](j, k) = mo(i * LIBMESH_DIM + j) * no(i * LIBMESH_DIM + k);
	
  // Output slip plane normals
  // slip_plane_normals includes the normals of all the slip systems
  for (unsigned int i = 0; i < _variable_size; ++i) {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j) {  
      slip_plane_normals[i * LIBMESH_DIM + j] = no(i * LIBMESH_DIM + j);		  
    }
  }  
  
}
