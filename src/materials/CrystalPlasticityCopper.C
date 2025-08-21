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
  params.addParam<Real>("h",0.0,"Direct hardening coefficient for backstress");
  params.addParam<Real>("h_D",0.0,"Dynamic recovery coefficient for backstress");
  return params;
}

CrystalPlasticityCopper::CrystalPlasticityCopper(
    const InputParameters & parameters)
  : CrystalPlasticityKalidindiUpdate(parameters),
    _gss_initial_std(getParam<Real>("gss_initial_std")),
    
    // Backstress parameters
    _h(getParam<Real>("h")),
    _h_D(getParam<Real>("h_D")),

    // Backstress variable
    _backstress(declareProperty<std::vector<Real>>("backstress")),
    _backstress_old(getMaterialPropertyOld<std::vector<Real>>("backstress"))
{
}

void
CrystalPlasticityCopper::initQpStatefulProperties()
{
  // This assigns _gss_initial to _slip_resistance
  CrystalPlasticityKalidindiUpdate::initQpStatefulProperties();
  
  Real rand_num;
  
  // Initialize the backstress size
  _backstress[_qp].resize(_number_slip_systems);
  
  for (const auto i : make_range(_number_slip_systems))
  {
    // Random number between 0 and 1
    rand_num = MooseRandom::rand();
    
    _slip_resistance[_qp][i] += _gss_initial_std * (2.0 * rand_num - 1.0);
    // TO DO: add check that slip resistance remains positive
    
    // Initialize backstress
    _backstress[_qp][i] = 0.0;
  }
}

void
CrystalPlasticityCopper::setInitialConstitutiveVariableValues()
{
  // Initialize state variables with the value at the previous time step
  _backstress[_qp] = _backstress_old[_qp];
  _previous_substep_backstress = _backstress_old[_qp];
}

void
CrystalPlasticityCopper::setSubstepConstitutiveVariableValues()
{
  // Inialize state variable of the next substep
  // with the value at the previous substep
  _backstress[_qp] = _previous_substep_backstress;
}
