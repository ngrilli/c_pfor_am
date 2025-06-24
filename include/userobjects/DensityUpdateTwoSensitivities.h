// Nicolò Grilli
// Università di Bristol
// 24 Giugno 2025

#pragma once

#include "ElementUserObject.h"
#include "MooseTypes.h"

/**
 * Element user object that performs SIMP optimization using a bisection algorithm using a volume
 * constraint. Two sensitivities are required.
 * Specifically developed for multi-objective optimization of a mechanical + a thermomechanical simulation.
 */
class DensityUpdateTwoSensitivities : public ElementUserObject
{
public:
  static InputParameters validParams();

  DensityUpdateTwoSensitivities(const InputParameters & parameters);

  virtual void initialize() override{};
  virtual void timestepSetup() override;
  virtual void execute() override;
  virtual void finalize() override{};
  virtual void threadJoin(const UserObject &) override{};

protected:
  /// The system mesh
  const MooseMesh & _mesh;
  /// The name of the pseudo-density variable
  const VariableName _design_density_name;
  /// The elasticity compliance sensitivity name
  const VariableName _density_sensitivity_name;
  /// The pseudo-density variable
  MooseWritableVariable * _design_density;
  /// The filtered density sensitivity variable
  const MooseWritableVariable * _density_sensitivity;
  /// The volume fraction to be enforced
  const Real _volume_fraction;
  /// The thermo-elasticity compliance sensitivity (second sensitivity) name
  /// and corresponding filtered sensitivity variable
  const VariableName _tm_sensitivity_name;
  const MooseWritableVariable * _tm_sensitivity;
  
private:
  struct ElementData
  {
    Real old_density;
    Real sensitivity; // first sensitivity
    Real tm_sensitivity; // second sensitivity
    Real volume;
    Real new_density;
    ElementData() = default;
    ElementData(Real dens, Real sens, Real tm_sens, Real vol, Real filt_dens)
      : old_density(dens), sensitivity(sens), tm_sensitivity(tm_sens), volume(vol), new_density(filt_dens)
    {
    }
  };

  /**
   * Gathers element data necessary to perform the bisection algorithm for optimization
   */
  void gatherElementData();

  /**
   * Performs the optimality criterion loop (bisection)
   */
  void performOptimCritLoop();

  Real computeUpdatedDensity(Real current_density, Real dc, Real tm_sensitivity, Real lmid);

  /// Total volume allowed for volume contraint
  Real _total_allowable_volume;

  /// Data structure to hold old density, 2 sensitivities, volume, current density.
  std::map<dof_id_type, ElementData> _elem_data_map;
  
  /// Lower bound for bisection algorithm
  const Real _lower_bound;
  /// Upper bound for bisection algorithm
  const Real _upper_bound;
};
