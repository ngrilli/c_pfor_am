// Marco Pelegatti
// Università di Udine
// Nicolò Grilli
// Università di Bristol
// 20 Dicembre 2025

#include "DislocationStructuresCyclicJump.h"
#include "libmesh/int_range.h"
#include "libmesh/utility.h"
#include <cmath>
#include "Function.h"

registerMooseObject("c_pfor_amApp", DislocationStructuresCyclicJump);

InputParameters
DislocationStructuresCyclicJump::validParams()
{
  InputParameters params = CrystalPlasticityCyclicDislocationStructures::validParams();
  params.addClassDescription("Dislocation based model for crystal plasticity "
                             "using the cyclic jump for acceleration. ");
  params.addCoupledVar("cyclic_jump", 0, "An auxiliary variable that indicates when cyclic jump update holds: "
                                         "0 = physics-based update; 1 = cyclic jump extrapolation");
  params.addParam<Real>("real_cycle_duration", 0, "How long do the real cycles last?");
  return params;
}

DislocationStructuresCyclicJump::DislocationStructuresCyclicJump(
    const InputParameters & parameters)
  : CrystalPlasticityCyclicDislocationStructures(parameters),
  _cyclic_jump(coupledValue("cyclic_jump")),
  _cyclic_jump_old(coupledValueOld("cyclic_jump")),
  _real_cycle_duration(getParam<Real>("real_cycle_duration")),
  
  // Extrapolated variables store the value at the beginning of the physical cycle
  // when the cyclic jump variable changes from 0 to 1 
  _extrapolated_slip_increment_c(declareProperty<std::vector<Real>>("slip_rate_c")),
  _extrapolated_slip_increment_w(declareProperty<std::vector<Real>>("slip_rate_w")),
  _extrapolated_slip_increment_PSB(declareProperty<std::vector<Real>>("slip_rate_PSB"))
  
  /// UTILE AVERE ANCHE LA OLD COSI' LA CONDIZIONE CHE OLD E NEW SONO DIVERSE
  /// TI DA LA CONDIZIONE PER CUI LE EXTRAPOLATED STATE VARIABLES VENGONO ASSEGNATE
  /// SICCOME TUTTE LE FUNZIONI IN QUESTA CLASSE VENGONO CHIAMATE PIU' VOLTE PER TIMESTEP
  /// NON C'E' UN MODO CHE SIA PIU' GIUSTO PER ASSEGNARE THE EXTRAPOLATED STATE VARIABLES
  /// CONSIGLIO DI METTERE UNA FUNZIONE DENTRO A calculateSlipRate
  /// FATTA IN MODO CHE CHIAMATA PIU' VOLTE PER TIMESTEP NON INTERFERISCA CON LE STATE VARIABLES NORMALI
{
}

void
DislocationStructuresCyclicJump::initQpStatefulProperties()
{
  CrystalPlasticityCyclicDislocationStructures::initQpStatefulProperties();
  
  // Initialize extrapolated slip increment size
  _extrapolated_slip_increment_c[_qp].resize(_number_slip_systems);
  _extrapolated_slip_increment_w[_qp].resize(_number_slip_systems);
  _extrapolated_slip_increment_PSB[_qp].resize(_number_slip_systems);
 
  for (const auto i : make_range(_number_slip_systems))
  {
    _extrapolated_slip_increment_c[_qp][i] = 0.0;
    _extrapolated_slip_increment_w[_qp][i] = 0.0;
    _extrapolated_slip_increment_PSB[_qp][i] = 0.0;
  }
}

/// setInitialConstitutiveVariableValues NON DOVREBBE ESSERE NECESSARIA
/// PER LE PROJECTED STATE VARIABLES
/// SOLO LE STATE VARIABLES NORMALI SEGUONO LA PROCEDURA DI SUBSTEPPING COME PRIMA
/// ASSUMENDO UPDATE BASATO SUL RATE CALCOLATO CON PROJECTED STATE VARIABLES
/// IN CASO DI CYCLIC JUMP
/// STESSA COSA PER setSubstepConstitutiveVariableValues
/// E PER updateSubstepConstitutiveVariableValues
/// E PER cacheStateVariablesBeforeUpdate

/// calculateSlipResistance SEMBRA CALCOLI SOLO VARIABILI DIPENDENTI
/// E SI PUO USARE LA STESSA DURANTE CYCLIC JUMP

// Slip resistance can be calculated from dislocation density here only
// because it is the first method in which it is used,
// while calculateConstitutiveSlipDerivative is called afterwards
bool
DislocationStructuresCyclicJump::calculateSlipRate()
{
  if (_cyclic_jump[_qp] > 0.5) {
	  
    /// QUI UPDATE BASATO SU ESTRAPOLAZIONE
    /// EXTRAPOLATED_RATE PUO ESSERE CALCOLATO QUI DIRETTAMENTE IN BASE
    /// ALL'INTERVALLO DI TEMPO PREDETERMINATO O CALCOLATO A RUNTIME 
    /// SUGGERISCO IMPORTARE FUNCTIONS PER L'INTERVALLO DI TEMPO
    /// SICCOME TI PERMETTE POI DI RENDERLO ADATTIVO
    /// _slip_increment_c[_qp][i] = EXTRAPOLATED_RATE_C
    /// _slip_increment_w[_qp][i] = EXTRAPOLATED_RATE_W
    
    for (const auto i : make_range(_number_slip_systems))
    {	  
	  // CHECK + NEED TO DEFINE OLD 
      //_slip_increment_c[_qp][i] = _slip_increment_c_old[_qp][i] + (_slip_increment_c_old[_qp][i] - _extrapolated_slip_increment_c[_qp][i]) * ((_t_old - _time_begin_real_cycles) / _real_cycle_duration);
    }
    
    /// LA STESSA COSA VA FATTA IN TUTTE LE FUNZIONI
    /// DOVE SI CALCOLANO GLI *_INCREMENT_*
    /// A QUEL PUNTO LA updateStateVariables
    /// E' GIA' CORRETTA E NON VA MODIFICATA
	  
  } else {
	  
    // condition for which the actual cycles start = cyclic jump variable changes from 1 to 0 
    // at this points, store the value that is required later for extrapolation
    // when the cyclic jump variable will change from 0 to 1
    if (_cyclic_jump_old[_qp] > 0.5) {
	
      for (const auto i : make_range(_number_slip_systems))
      {	
        _extrapolated_slip_increment_c[_qp][i] = _slip_increment_c[_qp][i];
        _extrapolated_slip_increment_w[_qp][i] = _slip_increment_w[_qp][i];
        _extrapolated_slip_increment_PSB[_qp][i] = _slip_increment_PSB[_qp][i];
        
        // time at which the last memorization has taken place
        _time_begin_real_cycles = _t;
      }
	}
	  
    CrystalPlasticityCyclicDislocationStructures::calculateSlipRate();
  }

  return true;
}
