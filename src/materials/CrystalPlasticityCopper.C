// Nicolò Grilli
// Università di Bristol
// 2 Agosto 2025

#include "CrystalPlasticityCopper.h"
#include "libmesh/int_range.h"
#include "MooseRandom.h"

registerMooseObject("c_pfor_amApp", CrystalPlasticityCopper);

InputParameters
CrystalPlasticityCopper::validParams()
{
  InputParameters params = CrystalPlasticityKalidindiUpdate::validParams();
  params.addClassDescription("Phenomenological model for copper with spatially random gss_initial.");
  params.addParam<Real>("gss_initial_std", 0.0, "Standard deviation of the spatial distribution of the " 
                                                "initial lattice friction strength of the material.");
  return params;
}

CrystalPlasticityCopper::CrystalPlasticityCopper(
    const InputParameters & parameters)
  : CrystalPlasticityKalidindiUpdate(parameters),
    _gss_initial_std(getParam<Real>("gss_initial_std"))
{
}

void
CrystalPlasticityCopper::initQpStatefulProperties()
{
  // This assigns _gss_initial to _slip_resistance
  CrystalPlasticityKalidindiUpdate::initQpStatefulProperties();
  
  Real rand_num;
  
  for (const auto i : make_range(_number_slip_systems))
  {
    // Random number between 0 and 1
    rand_num = MooseRandom::rand();
    
    _slip_resistance[_qp][i] += _gss_initial_std * (2.0 * rand_num - 1.0);
    // TO DO: add check that slip resistance remains positive
  }
}
