// Nicolò Grilli
// University of Bristol
// 30 Aprile 2023

#include "UMATStressR5CreepDamage.h"

registerMooseObject("TensorMechanicsApp", UMATStressR5CreepDamage);

InputParameters
UMATStressR5CreepDamage::validParams()
{
  InputParameters params = UMATStressDamage::validParams();
  params.addClassDescription("Coupling material to use Abaqus UMAT models in MOOSE "
                             "and coupling with phase field damage model. "      
                             "This is a variant of UMATStressDamage in which the fracture energy "
                             "is degraded based on the R5 creep damage criterion, which is based on "
                             "ductility exhaustion theory."
	                         "The damage model is reported in: "
	                         "M. W. Spindler, "
	                         "The prediction of creep damage in type 347 weld metal. "
	                         "Part I: the determination of material properties from creep and tensile tests "
	                         "International Journal of Pressure Vessels and Piping "
	                         "Volume 82, Issue 3, March 2005, Pages 175-184 ");
  params.addParam<Real>("residual_creep_degradation", 1e-3, "Minimum residual creep degradation.");
  params.addParam<MaterialPropertyName>("creep_degradation_name", "creep_degradation", "Creep degradation function.");
  return params;
}

UMATStressR5CreepDamage::UMATStressR5CreepDamage(const InputParameters & parameters)
  : UMATStressDamage(parameters),
  _residual_creep_degradation(getParam<Real>("residual_creep_degradation")),
  _creep_degradation_old(getMaterialPropertyOld<Real>("creep_degradation")),
  _f_ep_c(declareProperty<Real>("f_ep_c"))
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

// compute history variable and assign to _E
// which is used by the fracture model for damage growth
// Damage grows only because of the positive part of the elastic energy F_pos
// Creep damage degradation function is added, which corresponds to a decrease of Gc
// in the damage model
void
UMATStressR5CreepDamage::computeHistoryVariable(Real & F_pos, Real & F_neg)
{
  // Assign history variable
  Real hist_variable = _H_old[_qp];
  
  // Fpos as modified during creep damage
  // It is effectively equivalent to a decrease in Gc
  Real creep_F_pos;
  
  // Pass the creep damage function from UMAT to MOOSE material property
  _f_ep_c[_qp] = _aqSTATEV[9];
  
  // Calculation of _creep_degradation_old is done in the input file
  // using a ParsedMaterial as a function of _f_ep_c
  
  if (_creep_degradation_old[_qp] > _residual_creep_degradation)
    creep_F_pos = F_pos / _creep_degradation_old[_qp];
  else 
    creep_F_pos = F_pos / _residual_creep_degradation;
  
  // _use_snes_vi_solver option not implemented

  if (creep_F_pos > _H_old[_qp])
    _H[_qp] = creep_F_pos;
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
