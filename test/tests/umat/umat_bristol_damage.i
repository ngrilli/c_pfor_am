# Testing the UMAT Interface for crystal plasticity using finite strain deformation.
# The Bristol_CP code has a path inside the fortran code BRISTOL.f
# that must be changed to point to this folder
# would be better to read the path from inputs.dat in the future

# coupling with phase field damage model
# damage is constant here

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 2.0
    zmin = 0.0
    zmax = 2.0
	nx = 2
	ny = 2
	nz = 2
  []
  [origin_modifier]
    type = BoundingBoxNodeSetGenerator
    input = gen
    new_boundary = origin
    top_right = '0.01 0.01 0.01'
    bottom_left = '-0.01 -0.01 -0.01'
  []  
  [xp_point_modifier]
    type = BoundingBoxNodeSetGenerator
    input = origin_modifier
    new_boundary = xp_point
    top_right = '2.01 0.01 0.01'
    bottom_left = '1.99 -0.01 -0.01'
  []  
  [zp_point_modifier]
    type = BoundingBoxNodeSetGenerator
    input = xp_point_modifier
    new_boundary = zp_point
    top_right = '0.01 0.01 1.01'
    bottom_left = '-0.01 -0.01 0.99'
  []  
[]

# add damage variable
# it is constant here
[Variables]
  [./c]
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = '0.5'
    [../]
  [../]
[]

[Functions]
  [top_pull]
    type = ParsedFunction
    value = 't/1000'
  []
  [lateral_pressure]
    type = ParsedFunction
	value = '2.0'
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    add_variables = true
    strain = FINITE
  []
[]

# constant damage
[Kernels]
  [./dcdt]
    type = TimeDerivative
    variable = c
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
  [./statev_1]
    order = CONSTANT
    family = MONOMIAL  
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
  [./statev_1]
    type = MaterialStdVectorAux
    variable = statev_1
    property = state_var
    index = 0
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
  [Pressure]
    [bc_pressure_front]
      boundary = front
      function = lateral_pressure
    []
    [bc_pressure_right]
      boundary = right
      function = lateral_pressure
    []
  []
[]

# UMAT interface, there are no STATEV used
# because fortran modules are adopted
# modules need to be compiled independently
# with gfortran before make
# command is for instance:
# gfortran -c -free calculations.f
[Materials]
  [umat]
    type = UMATStressDamage
    constant_properties = '0'
    plugin = '../../plugins/Bristol_CP/BRISTOL'
    num_state_vars = 9
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
    constant_expressions = '1.0e-6'
    derivative_order = 2
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '2.0 0.5 1e-6'
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
  
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  
  line_search = 'none'

  nl_max_its = 100
  l_max_its = 10
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  
  start_time = 0.0
  end_time = 0.01 # run until 18.0 to see the elasto-plastic stress strain curve
  
  # if the number of non-linear iterations is in the interval
  # [optimal_iterations-iteration_window ; optimal_iterations+iteration_window]
  # then time step is preserved, otherwise it is cut
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25
    iteration_window = 15
    growth_factor = 1.0001
  	reject_large_step = True
    cutback_factor = 0.5
    timestep_limiting_postprocessor = matl_ts_min
    dt = 0.01
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

