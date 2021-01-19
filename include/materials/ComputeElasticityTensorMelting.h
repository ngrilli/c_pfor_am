// Nicolo Grilli
// National University of Singapore
// 13 Novembre 2020

#pragma once

#include "ComputeElasticityTensorCPGrain.h"
#include "LaserTempReadFile.h"

/**
 * ComputeElasticityTensorMelting defines an elasticity tensor material object for crystal plasticity.
 * It is based on the user object GrainPropertyReadFile
 * that allows to assign grains from the euler angles input file
 * to the physical volumes in GMSH
 * Melting is considered: stiffness is degraded when the
 * temperature increases above melting or below gas temperature.
 */
class ComputeElasticityTensorMelting : public ComputeElasticityTensorCPGrain
{
public:
  static InputParameters validParams();

  ComputeElasticityTensorMelting(const InputParameters & parameters);

protected:

  virtual void computeQpElasticityTensor() override;
  
  virtual void checkPhase();
  
  virtual void melting();
  
  const Real _melting_temperature_high;
  const Real _melting_temperature_low;
  const Real _gas_temperature_high;
  const Real _gas_temperature_low;
  
  /// Residual stiffness of gas and molten pool (percent)
  const Real _residual_stiffness;

  /// The LaserTempReadFile GeneralUserObject to read element specific temperature values from file
  const LaserTempReadFile * const _temperature_read_user_object;
  
  /// Time interval between two temperature data field
  const Real _temperature_time_step;
  
  /// Using the element activation feature
  /// Elements must be activated with degraded stiffness
  bool _activate_elems;  

  /// Stiffness tensor modified by melting
  /// to model laser scanning
  RankFourTensor _Melt_Cijkl;
  
  /// Flags to indicate the phase at the current 
  /// and next temperature time step 
  unsigned int _isSolid;
  unsigned int _isLiquid;
  unsigned int _isGas;
  unsigned int _isSolidNext;
  unsigned int _isLiquidNext;
  unsigned int _isGasNext;
  
  /// Flags to indicate the phase at the previous temperature time step 
  /// to understand if element was activated at this temperature time step
  unsigned int _isSolidPrevious;  
  
  // Temperature values at the current and 
  // next temperature time step
  Real _TempValue;
  Real _TempValueNext;
  
  // Temperature value at the previous temperature time step
  Real _TempValuePrevious;
  
  // temperature time step to be used from the CFD simulations
  unsigned int _temperature_step;
  
  // fraction of temperature time step completed, between 0 and 1
  Real _FracTimeStep;
  
};
