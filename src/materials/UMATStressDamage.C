// Nicolò Grilli
// University of Bristol
// 26 Gennaio 2022

#include "UMATStressDamage.h"
#include "Factory.h"
#include "MooseMesh.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "libmesh/int_range.h"
#include "petscblaslapack.h"
#include "MooseException.h"
#include "MathUtils.h"
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
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>(
      "E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<Real>("bulk_modulus_ref",0.0,"reference bulk modulus for vol/non-vol decomposition");
  return params;
}

#ifndef METHOD
#error "The METHOD preprocessor symbol must be supplied by the build system."
#endif

UMATStressDamage::UMATStressDamage(const InputParameters & parameters)
  : AbaqusUMATStress(parameters),
    _c(coupledValue("c")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _H(declareProperty<Real>("hist")), // History variable to avoid damage decrease 
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _E(declareProperty<Real>(getParam<MaterialPropertyName>("E_name"))),
    _dEdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"),
                                          getVar("c", 0)->name())),
    _d2Ed2c(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _dstress_dc(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())), 
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),		
	_D(getMaterialProperty<Real>("D_name")),
    _dDdc(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name())),
    _d2Dd2c(getMaterialPropertyDerivative<Real>(
        "D_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _bulk_modulus_ref(getParam<Real>("bulk_modulus_ref")) // reference bulk modulus for vol/non-vol decomposition
{
  if (!_use_one_based_indexing)
    mooseDeprecated(
        "AbaqusUMATStress has transitioned to 1-based indexing in the element (NOEL) and "
        "integration point (NPT) numbers to ensure maximum compatibility with legacy UMAT files. "
        "Please ensure that any new UMAT plugins using these quantities are using the correct "
        "indexing. 0-based indexing will be deprecated soon.");	
	
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

  Real F_pos, F_neg; // tensile and compressive part of the elastic strain energy

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
  // The order of the strain components is adapted to the one in Abaqus
  static const std::array<Real, 6> strain_factor{{1, 1, 1, 2, 2, 2}};
  static const std::array<std::pair<unsigned int, unsigned int>, 6> component{
      {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}}};

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
  
  // integration point number
  _aqNPT = _qp + (_use_one_based_indexing ? 1 : 0);
  
  // damage phase field is passed to STATEV(1) 
  _aqSTATEV[0] = _c[_qp];

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
        &_aqNPT,
        &_aqLAYER,
        &_aqKSPT,
        &_aqKSTEP,
        &_aqKINC);

  // Update state variables
  for (int i = 0; i < _aqNSTATV; ++i)
    _state_var[_qp][i] = _aqSTATEV[i];

  // all the quantities necessary for
  // the phase field damage model are defined
  assignFreeEnergy(F_pos, F_neg);
  
  // Compute history variables necessary for
  // the phase field damage model
  computeHistoryVariable(F_pos, F_neg);

  // Here, we apply UMAT convention: Always multiply _dt by PNEWDT to determine the material time
  // step MOOSE time stepper will choose the most limiting of all material time step increments
  // provided
  _material_timestep[_qp] = _aqPNEWDT * _dt;
	  
  // Get new stress tensor - UMAT should update stress
  _stress[_qp] = RankTwoTensor(
      _aqSTRESS[0], _aqSTRESS[1], _aqSTRESS[2], _aqSTRESS[5], _aqSTRESS[4], _aqSTRESS[3]);

  // Rotate the stress state to the current configuration
  _stress[_qp].rotate(_rotation_increment[_qp]);
  
  // Build Jacobian matrix from UMAT's Voigt non-standard order to fourth order tensor.
  
  const unsigned int N = LIBMESH_DIM;
  const unsigned int ntens = N * (N + 1) / 2;
  const int nskip = N - 1;

  for (auto i : make_range(N))
    for (auto j : make_range(N))
      for (auto k : make_range(N))
        for (auto l : make_range(N))
        {
          if (i == j)
            _jacobian_mult[_qp](i, j, k, l) =
                k == l ? _aqDDSDDE[i * ntens + k] : _aqDDSDDE[i * ntens + k + nskip + l];
          else
            // i!=j
            _jacobian_mult[_qp](i, j, k, l) =
                k == l ? _aqDDSDDE[(nskip + i + j) * ntens + k]
                       : _aqDDSDDE[(nskip + i + j) * ntens + k + nskip + l];
        }

}

// Free energy components and their derivatives
// calculated by the UMAT
// are assigned to the corresponding MaterialProperty
void
UMATStressDamage::assignFreeEnergy(Real & F_pos, Real & F_neg)
{
  // Positive part of the second Piola-Kirchhoff stress
  RankTwoTensor pk2_pos;	
	
  // Assign positive and negative parts of the free energy
  // Equations 13 and 14 in Grilli, Koslowski, 2019
  F_pos = _aqSTATEV[1];
  F_neg = _aqSTATEV[2];
  
  // Assign positive part of the second Piola-Kirchhoff stress
  // it follows Abaqus components order convention for 6-vectors
  // Therefore fourth and sixth components must be swapped 
  pk2_pos = RankTwoTensor(
      _aqSTATEV[3], _aqSTATEV[4], _aqSTATEV[5], _aqSTATEV[8], _aqSTATEV[7], _aqSTATEV[6]);
  
  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = pk2_pos * _dDdc[_qp];
  
  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history variable
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = pk2_pos * _dDdc[_qp];
  
}

// compute history variable and assign to _E
// which is used by the fracture model for damage growth
// Damage grows only because of the positive part of the elastic energy F_pos
void
UMATStressDamage::computeHistoryVariable(Real & F_pos, Real & F_neg)
{
  // Assign history variable
  Real hist_variable = _H_old[_qp];
  
  // _use_snes_vi_solver option not implemented

  if (F_pos > _H_old[_qp])
    _H[_qp] = F_pos;
  else
    _H[_qp] = _H_old[_qp];

  if (_use_current_hist)
    hist_variable = _H[_qp];

  // _barrier not implemented

  // Elastic free energy density and derivatives
  _E[_qp] = hist_variable * _D[_qp] + F_neg;
  _dEdc[_qp] = hist_variable * _dDdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp];

}
