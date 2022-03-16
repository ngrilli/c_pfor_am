# Testing the UMAT Interface for crystal plasticity using finite strain deformation.
# coupling with phase field damage model
# damage evolution is implemented

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 1.0
	nx = 1
	ny = 1
	nz = 1
  []
[]

[Functions]
  [top_pull]
    type = ParsedFunction
    value = 't/1000'
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    add_variables = true
    strain = FINITE
  []
[]

[Modules/PhaseField/Nonconserved/c]
  free_energy = F
  kappa = kappa_op
  mobility = L
[] 

# off diagonal Jacobian kernels
[Kernels]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
    use_displaced_mesh = true
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
    use_displaced_mesh = true
  [../]
  [./solid_z]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_z
    component = 2
    c = c
    use_displaced_mesh = true
  [../]
  [./off_disp]
    type = AllenCahnElasticEnergyOffDiag
    variable = c
    displacements = 'disp_x disp_y disp_z'
    mob_name = L
    use_displaced_mesh = true
  [../]
[]

[AuxVariables]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./stress_yx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./stress_zx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./stress_zy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  # this state variable is a constant monomial
  # copy of the damage variable
  [./state_damage]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  
  # positive part of the free energy
  [./Fpos]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  
  # positive part of the free energy
  [./Fneg]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  
  # components of the tensile part
  # of the second Piola-Kirchhoff stress
  [./pk2_pos_xx]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  [./pk2_pos_yy]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  [./pk2_pos_zz]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  [./pk2_pos_xy]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  [./pk2_pos_xz]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  [./pk2_pos_yz]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  
  [./bounds_dummy]
  [../]

[]

# bounds used to impose that damage can only grow
[Bounds]
  [./irreversibility]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = lower
  [../]
  [./upper]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = upper
    bound_value = 1.0
  [../]
[]

[AuxKernels] 

  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]  
  
  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  
  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  [../]
  
  [./stress_xz]
    type = RankTwoAux
    variable = stress_xz
    rank_two_tensor = stress
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  [../]

  [./stress_yx]
    type = RankTwoAux
    variable = stress_yx
    rank_two_tensor = stress
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  [../]  

  [./stress_yz]
    type = RankTwoAux
    variable = stress_yz
    rank_two_tensor = stress
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  [../]  

  [./stress_zx]
    type = RankTwoAux
    variable = stress_zx
    rank_two_tensor = stress
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  [../]  

  [./stress_zy]
    type = RankTwoAux
    variable = stress_zy
    rank_two_tensor = stress
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  [../]
  
  # this state variable in a constant mononial
  # copy of the damage variable
  [./state_damage]
    type = MaterialStdVectorAux
    variable = state_damage
    property = state_var
    index = 0
    execute_on = timestep_end
  [../]
  
  [./Fpos]
    type = MaterialStdVectorAux
    variable = Fpos
    property = state_var
    index = 1
    execute_on = timestep_end  
  [../]
  
  [./Fneg]
    type = MaterialStdVectorAux
    variable = Fneg
    property = state_var
    index = 2
    execute_on = timestep_end  
  [../]
  
  [./pk2_pos_xx]
    type = MaterialStdVectorAux
    variable = pk2_pos_xx
    property = state_var
    index = 3
    execute_on = timestep_end    
  [../]
  
  [./pk2_pos_yy]
    type = MaterialStdVectorAux
    variable = pk2_pos_yy
    property = state_var
    index = 4
    execute_on = timestep_end    
  [../]
  
  [./pk2_pos_zz]
    type = MaterialStdVectorAux
    variable = pk2_pos_zz
    property = state_var
    index = 5
    execute_on = timestep_end    
  [../]
  
  [./pk2_pos_xy]
    type = MaterialStdVectorAux
    variable = pk2_pos_xy
    property = state_var
    index = 6
    execute_on = timestep_end    
  [../]
  
  [./pk2_pos_xz]
    type = MaterialStdVectorAux
    variable = pk2_pos_xz
    property = state_var
    index = 7
    execute_on = timestep_end    
  [../]
  
  [./pk2_pos_yz]
    type = MaterialStdVectorAux
    variable = pk2_pos_yz
    property = state_var
    index = 8
    execute_on = timestep_end    
  [../]


[]

# pure tension boundary conditions
[BCs]
  [y_pull_function]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = top_pull
  []
  [x_bot]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  []
  [y_bot]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []
  [z_bot]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  []
[]

[Materials]
  [umat]
    type = UMATStressDamage
    constant_properties = '0'
    plugin = '../../plugins/Bristol_CP/BRISTOL'
    num_state_vars = 72
	c = c
	use_one_based_indexing = true
    E_name = 'elastic_energy'
    D_name = 'degradation'
    use_current_history_variable = true
  []
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '1.0e-5'
    derivative_order = 2
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '2.0 4.0 5e-6'
  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l'
  [../]
  [./local_fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'c'
    material_property_names = 'gc_prop l'
    function = 'c^2 * gc_prop / 2 / l'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    args = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    f_name = F
  [../]
[]

[UserObjects]

  # this runs the UExternalDB subroutine once at the first time step
  [uexternaldb]
    type = AbaqusUExternalDB
    plugin = '../../plugins/Bristol_CP/BRISTOL'
    execute_on = 'INITIAL'
  []
  
  # this is used to find unconverged time steps depending on
  # the UMAT output variable PNEWDT
  [time_step_size]
    type = TimestepSize
    execute_on = 'INITIAL LINEAR'
  []
  [terminator_umat]
    type = Terminator
    expression = 'time_step_size > matl_ts_min'
    fail_mode = SOFT
    execute_on = 'FINAL'
  []
[]

[Postprocessors]
  [matl_ts_min]
    type = MaterialTimeStepPostprocessor
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'
  
  line_search = 'none'

  nl_max_its = 20
  l_max_its = 20
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  
  start_time = 0.0
  end_time = 0.001 # run until 10.0 to see the elasto-plastic-damage stress strain curve
  
  # if the number of non-linear iterations is in the interval
  # [optimal_iterations-iteration_window ; optimal_iterations+iteration_window]
  # then time step is preserved, otherwise it is cut
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25
    iteration_window = 15
    growth_factor = 1.01
  	reject_large_step = True
    cutback_factor = 0.5
    timestep_limiting_postprocessor = matl_ts_min
    dt = 0.001
  []
  
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Outputs]
  exodus = true
  interval = 1
[]

