// Nicolò Grilli
// University of Bristol
// 19 Gennaio 2023

#include "HEPFlowRatePowerLawJ2.h"

registerMooseObject("TensorMechanicsApp", HEPFlowRatePowerLawJ2);

InputParameters
HEPFlowRatePowerLawJ2::validParams()
{
  InputParameters params = HEVPFlowRatePowerLawJ2::validParams();
  params.addClassDescription(
      "User object to evaluate power law flow rate and flow direction based on J2. "
      "It is similar to HEVPFlowRatePowerLawJ2 but "
      "plastic and not visco-plastic, therefore it gives zero "
      "plastic strain rate below the yield point. ");
  params.addParam<Real>("threshold_stress_ratio", 0.0, "Ratio between equivalent stress and yield strength "
                                                       "at which plastic slip starts. ");
  return params;
}

HEPFlowRatePowerLawJ2::HEPFlowRatePowerLawJ2(const InputParameters & parameters)
  : HEVPFlowRatePowerLawJ2(parameters),
    _threshold_stress_ratio(getParam<Real>("threshold_stress_ratio"))
{
}

bool
HEPFlowRatePowerLawJ2::computeValue(unsigned int qp, Real & val) const
{
  RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
  Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
  
  Real stress_ratio = eqv_stress / _strength[qp];
  
  if (stress_ratio < 1.0) { // below yield point
	  
    val = 0.0;
  
  } else { // above yield point
	  
    val = std::pow(stress_ratio, _flow_rate_exponent) * _ref_flow_rate;
	  
  }
  
  if (val > _flow_rate_tol)
  {
#ifdef DEBUG
    mooseWarning(
        "Flow rate greater than ", _flow_rate_tol, " ", val, " ", eqv_stress, " ", _strength[qp]);
#endif
    return false;
  }
  return true;
}

bool
HEPFlowRatePowerLawJ2::computeDerivative(unsigned int qp,
                                          const std::string & coupled_var_name,
                                          Real & val) const
{
  val = 0.0;

  if (_strength_prop_name == coupled_var_name)
  {
    RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
    Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
    
    Real stress_ratio = eqv_stress / _strength[qp];
    
    if (stress_ratio < 1.0) { // below yield point
		
      val = 0.0;
		
    } else { // above yield point
	
      val = -_ref_flow_rate * _flow_rate_exponent *
            std::pow(eqv_stress / _strength[qp], _flow_rate_exponent) / _strength[qp];
	
    }
  }

  return true;
}

bool
HEPFlowRatePowerLawJ2::computeTensorDerivative(unsigned int qp,
                                                const std::string & coupled_var_name,
                                                RankTwoTensor & val) const
{
  val.zero();

  if (_pk2_prop_name == coupled_var_name)
  {
    RankTwoTensor pk2_dev = computePK2Deviatoric(_pk2[qp], _ce[qp]);
    Real eqv_stress = computeEqvStress(pk2_dev, _ce[qp]);
    
    Real dflowrate_dseqv;
    
    Real stress_ratio = eqv_stress / _strength[qp];
    
    if (stress_ratio < 1.0) { // below yield point
		
      dflowrate_dseqv = 0.0;
		
    } else { // above yield point
	
      dflowrate_dseqv = _ref_flow_rate * _flow_rate_exponent *
                        std::pow(eqv_stress / _strength[qp], _flow_rate_exponent - 1.0) /
                        _strength[qp];
	
    }    

    RankTwoTensor tau = pk2_dev * _ce[qp];
    RankTwoTensor dseqv_dpk2dev;

    dseqv_dpk2dev.zero();
    if (eqv_stress > 0.0)
      dseqv_dpk2dev = 1.5 / eqv_stress * tau * _ce[qp];

    RankTwoTensor ce_inv = _ce[qp].inverse();

    RankFourTensor dpk2dev_dpk2;
    for (const auto i : make_range(Moose::dim))
      for (const auto j : make_range(Moose::dim))
        for (const auto k : make_range(Moose::dim))
          for (const auto l : make_range(Moose::dim))
          {
            dpk2dev_dpk2(i, j, k, l) = 0.0;
            if (i == k && j == l)
              dpk2dev_dpk2(i, j, k, l) = 1.0;
            dpk2dev_dpk2(i, j, k, l) -= ce_inv(i, j) * _ce[qp](k, l) / 3.0;
          }
    val = dflowrate_dseqv * dpk2dev_dpk2.transposeMajor() * dseqv_dpk2dev;
  }
  return true;
}

