// Nicolò Grilli
// University of Bristol
// 17 Marzo 2022

#include "CrystalPlasticityPosNegStateVariable.h"

#include <fstream>

registerMooseObject("TensorMechanicsApp", CrystalPlasticityPosNegStateVariable);

InputParameters
CrystalPlasticityPosNegStateVariable::validParams()
{
  InputParameters params = CrystalPlasticityStateVariable::validParams();
  params.addClassDescription(
      "Crystal plasticity state variable class. It can handle state variables "
      "that can become positive or negative during time. "	  
	  "Override the virtual functions in your class");
  return params;
}

CrystalPlasticityPosNegStateVariable::CrystalPlasticityPosNegStateVariable(const InputParameters & parameters)
  : CrystalPlasticityStateVariable(parameters)
{
  if (_scale_factor.size() != _num_mat_state_var_evol_rate_comps)
    mooseError("CrystalPlasticityPosNegStateVariable: Scale factor should be have the same size of "
               "evolution rate components.");

  _mat_prop_state_var_evol_rate_comps.resize(_num_mat_state_var_evol_rate_comps);

  for (unsigned int i = 0; i < _num_mat_state_var_evol_rate_comps; ++i)
    _mat_prop_state_var_evol_rate_comps[i] = &getMaterialProperty<std::vector<Real>>(
        parameters.get<std::vector<std::string>>("uo_state_var_evol_rate_comp_name")[i]);
}

// Update state variable, allow variables to become negative
// without limiting the lower value to _zero
bool
CrystalPlasticityPosNegStateVariable::updateStateVariable(unsigned int qp,
                                                    Real dt,
                                                    std::vector<Real> & val,
                                                    std::vector<Real> & val_old) const
{
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    val[i] = 0.0;
    for (unsigned int j = 0; j < _num_mat_state_var_evol_rate_comps; j++)
      val[i] += (*_mat_prop_state_var_evol_rate_comps[j])[qp][i] * dt * _scale_factor[j];
  }

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    val[i] = val_old[i] + val[i];
  }
  
  return true;
}
