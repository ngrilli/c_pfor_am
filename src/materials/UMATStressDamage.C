// Nicolò Grilli
// University of Bristol
// 19 Gennaio 2022

#include "UMATStressDamage.h"
#include "Factory.h"
#include "MooseMesh.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "libmesh/int_range.h"
#include <string.h>
#include <algorithm>

#define QUOTE(macro) stringifyName(macro)

registerMooseObject("TensorMechanicsApp", UMATStressDamage);

InputParameters
UMATStressDamage::validParams()
{
  InputParameters params = AbaqusUMATStress::validParams();
  params.addClassDescription("Coupling material to use Abaqus UMAT models in MOOSE "
                             "and coupling with phase field damage model");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  return params;
}

#ifndef METHOD
#error "The METHOD preprocessor symbol must be supplied by the build system."
#endif

UMATStressDamage::UMATStressDamage(const InputParameters & parameters)
  : AbaqusUMATStress(parameters),
    _c(coupledValue("c"))
{
  // get material properties
  for (std::size_t i = 0; i < _number_external_properties; ++i)
  {
    _external_properties[i] = &getMaterialProperty<Real>(_external_property_names[i]);
    _external_properties_old[i] = &getMaterialPropertyOld<Real>(_external_property_names[i]);
  }

  // Read mesh dimension and size UMAT arrays (we always size for full 3D)
  _aqNTENS = 6; // Size of the stress or strain component array (NDI+NSHR)
  _aqNSHR = 3;  // Number of engineering shear stress components
  _aqNDI = 3;   // Number of direct stress components (always 3)

  _aqDDSDDT.resize(_aqNTENS);
  _aqDRPLDE.resize(_aqNTENS);
  _aqSTRAN.resize(_aqNTENS);
  _aqDFGRD0.resize(9);
  _aqDFGRD1.resize(9);
  _aqDROT.resize(9);
  _aqSTRESS.resize(_aqNTENS);
  _aqDDSDDE.resize(_aqNTENS * _aqNTENS);
  _aqDSTRAN.resize(_aqNTENS);
  _aqPREDEF.resize(_number_external_fields + _number_external_properties);
  _aqDPRED.resize(_number_external_fields + _number_external_properties);
}


// initQpStatefulProperties needs to be checked
// because moose must be able to initialize
// damage variable without the material object affecting
// the IC object


void
UMATStressDamage::computeQpStress()
{
  const Real * myDFGRD0 = &(_Fbar_old[_qp](0, 0));
  const Real * myDFGRD1 = &(_Fbar[_qp](0, 0));
  const Real * myDROT = &(_rotation_increment[_qp](0, 0));

  // copy because UMAT does not guarantee constness
  for (unsigned int i = 0; i < 9; ++i)
  {
    _aqDFGRD0[i] = myDFGRD0[i];
    _aqDFGRD1[i] = myDFGRD1[i];
    _aqDROT[i] = myDROT[i];
  }

  // Recover "old" state variables
  for (int i = 0; i < _aqNSTATV; ++i)
    _aqSTATEV[i] = _state_var_old[_qp][i];

  // Pass through updated stress, total strain, and strain increment arrays
  static const std::array<Real, 6> strain_factor{{1, 1, 1, 2, 2, 2}};
  static const std::array<std::pair<unsigned int, unsigned int>, 6> component{
      {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}}};

  for (int i = 0; i < _aqNTENS; ++i)
  {
    _aqSTRESS[i] = _stress_old[_qp](component[i].first, component[i].second);
    _aqSTRAN[i] =
        _total_strain_old[_qp](component[i].first, component[i].second) * strain_factor[i];
    _aqDSTRAN[i] =
        _strain_increment[_qp](component[i].first, component[i].second) * strain_factor[i];
  }

  // current coordinates
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    _aqCOORDS[i] = _q_point[_qp](i);

  // zero out Jacobian contribution
  for (int i = 0; i < _aqNTENS * _aqNTENS; ++i)
    _aqDDSDDE[i] = 0.0;

  // Set PNEWDT initially to a large value
  _aqPNEWDT = std::numeric_limits<Real>::max();

  // Temperature
  _aqTEMP = _temperature_old[_qp];

  // Temperature increment
  _aqDTEMP = _temperature[_qp] - _temperature_old[_qp];

  for (const auto i : make_range(_number_external_fields))
  {
    // External field at this step
    _aqPREDEF[i] = (*_external_fields_old[i])[_qp];

    // External field increments
    _aqDPRED[i] = (*_external_fields[i])[_qp] - (*_external_fields_old[i])[_qp];
  }

  for (const auto i : make_range(_number_external_properties))
  {
    // External property at this step
    _aqPREDEF[i + _number_external_fields] = (*_external_properties_old[i])[_qp];

    // External property increments
    _aqDPRED[i + _number_external_fields] =
        (*_external_properties[i])[_qp] - (*_external_properties_old[i])[_qp];
  }

  // Layer number (not supported)
  _aqLAYER = -1;

  // Section point number within the layer (not supported)
  _aqKSPT = -1;

  // Increment number
  _aqKINC = _t_step;
  _aqKSTEP = 1;
  
  // damage needs to be passed to the state variable here

  // Connection to extern statement
  _umat(_aqSTRESS.data(),
        _aqSTATEV.data(),
        _aqDDSDDE.data(),
        &_elastic_strain_energy[_qp],
        &_plastic_dissipation[_qp],
        &_creep_dissipation[_qp],
        &_aqRPL,
        _aqDDSDDT.data(),
        _aqDRPLDE.data(),
        &_aqDRPLDT,
        _aqSTRAN.data(),
        _aqDSTRAN.data(),
        _aqTIME.data(),
        &_aqDTIME,
        &_aqTEMP,
        &_aqDTEMP,
        _aqPREDEF.data(),
        _aqDPRED.data(),
        _aqCMNAME,
        &_aqNDI,
        &_aqNSHR,
        &_aqNTENS,
        &_aqNSTATV,
        _aqPROPS.data(),
        &_aqNPROPS,
        _aqCOORDS.data(),
        _aqDROT.data(),
        &_aqPNEWDT,
        &_aqCELENT,
        _aqDFGRD0.data(),
        _aqDFGRD1.data(),
        &_aqNOEL,
        &_qp,
        &_aqLAYER,
        &_aqKSPT,
        &_aqKSTEP,
        &_aqKINC);

  // Update state variables
  for (int i = 0; i < _aqNSTATV; ++i)
    _state_var[_qp][i] = _aqSTATEV[i];

  // here you need to define all the quantities
  // necessary for the phase field model

  // Here, we apply UMAT convention: Always multiply _dt by PNEWDT to determine the material time
  // step MOOSE time stepper will choose the most limiting of all material time step increments
  // provided
  _material_timestep[_qp] = _aqPNEWDT * _dt;

  // Get new stress tensor - UMAT should update stress
  _stress[_qp] = RankTwoTensor(
      _aqSTRESS[0], _aqSTRESS[1], _aqSTRESS[2], _aqSTRESS[3], _aqSTRESS[4], _aqSTRESS[5]);

  // Rotate the stress state to the current configuration
  _stress[_qp].rotate(_rotation_increment[_qp]);

  // use DDSDDE as Jacobian mult
  _jacobian_mult[_qp].fillSymmetric21FromInputVector(std::array<Real, 21>{{
      _aqDDSDDE[0],  // C1111
      _aqDDSDDE[1],  // C1122
      _aqDDSDDE[2],  // C1133
      _aqDDSDDE[3],  // C1123
      _aqDDSDDE[4],  // C1113
      _aqDDSDDE[5],  // C1112
      _aqDDSDDE[7],  // C2222
      _aqDDSDDE[8],  // C2233
      _aqDDSDDE[9],  // C2223
      _aqDDSDDE[10], // C2213
      _aqDDSDDE[11], // C2212
      _aqDDSDDE[14], // C3333
      _aqDDSDDE[15], // C3323
      _aqDDSDDE[16], // C3313
      _aqDDSDDE[17], // C3312
      _aqDDSDDE[21], // C2323
      _aqDDSDDE[22], // C2313
      _aqDDSDDE[23], // C2312
      _aqDDSDDE[28], // C1313
      _aqDDSDDE[29], // C1312
      _aqDDSDDE[35]  // C1212
  }});
}
