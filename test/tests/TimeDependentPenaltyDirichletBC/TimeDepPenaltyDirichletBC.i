# test for TimeDependentPenaltyDirichletBC
# first a uniaxial tension with blocked lateral boundaries is applied
# from time 0 to 1, then from time 1 to 2
# the penalty function is decreased to zero
# and the lateral boundaries become free boundaries

[GlobalParams]
  displacements = 'ux uy uz'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
[]

[AuxVariables]
  [./fp_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    value = 'if(0.5*t,0.01,0.01*t)'
  [../]
  [./penalty_function]
    type = ParsedFunction
    value = '1e9*if(0.5*t,2-t,1)'
  [../]
  [./zero]
    type = ParsedFunction
    value = '0.0'
  [../]
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
  generate_output = stress_zz
[]

[AuxKernels]
  [./fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = plastic_deformation_gradient
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = total_lagrangian_strain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  [./gss1]
    type = MaterialStdVectorAux
    variable = gss
    property = slip_resistance
    index = 0
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./symmy]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0
  [../]
  [./symmx]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0
  [../]
  [./symmz]
    type = DirichletBC
    variable = uz
    boundary = back
    value = 0
  [../]
  [./tdisp]
    type = FunctionDirichletBC
    variable = uz
    boundary = front
    function = tdisp
  [../]
  [./penalty_right]
    type = TimeDependentPenaltyDirichletBC
    variable = ux
    boundary = right
    forcing_function = zero
    penalty_function = penalty_function
  [../]
  [./penalty_top]
    type = TimeDependentPenaltyDirichletBC
    variable = uy
    boundary = top
    forcing_function = zero
    penalty_function = penalty_function
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensorConstantRotationCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
  [../]
  [./stress]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    maximum_substep_iteration = 200
    use_line_search = true
    min_line_search_step_size = 0.01
  [../]
  [./trial_xtalpl]
    type = CrystalPlasticityKalidindiUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
    resistance_tol = 0.01
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.1
  solve_type = 'PJFNK'

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomerang
  nl_abs_tol = 1e-8
  nl_rel_step_tol = 1e-8
  dtmax = 0.1
  nl_rel_tol = 1e-8
  end_time = 2
  dtmin = 0.01
  nl_abs_step_tol = 1e-10
[]

[Outputs]
  exodus = true
  interval = 5
[]
