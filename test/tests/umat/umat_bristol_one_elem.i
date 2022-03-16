# Testing the UMAT Interface for crystal plasticity using finite strain deformation.
# The Bristol_CP code has a path inside the fortran code BRISTOL.f
# that must be changed to point to this folder
# would be better to read the path from inputs.dat in the future

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = -0.5
    xmax = 0.5
    ymin = -0.5
    ymax = 0.5
    zmin = -0.5
    zmax = 0.5
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
    boundary = front
    value = 0.0
  []
[]

# UMAT interface, there are no STATEV used
# because fortran modules are adopted
# modules need to be compiled independently
# with gfortran before make
# procedure is the following on first installation:
# convert BRISTOL.for into BRISTOL.f
# make sure all previous .mod, .o, .plugin files
# are removed, including the .mod files in the c_pfor_am folder
# gfortran -c -free calculations.f
# on each module, note BRISTOL.f is not a module
# make
# after code modification,
# remove .mod, .o, .plugin of the modified code
# remove BRISTOL.plugin
# remove .mod files in the c_pfor_am folder
# recompile modified module and make
[Materials]
  [umat]
    type = AbaqusUMATStress
    constant_properties = '0'
    plugin = '../../plugins/Bristol_CP/BRISTOL'
    num_state_vars = 72
	use_one_based_indexing = true
  []
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
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg          101'
  line_search = 'none'

  nl_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  
  start_time = 0.0
  end_time = 0.02 # run until 3.0 to see the elasto-plastic stress strain curve
  
  # if the number of non-linear iterations is in the interval
  # [optimal_iterations-iteration_window ; optimal_iterations+iteration_window]
  # then time step is preserved, otherwise it is cut
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25
    iteration_window = 20
    growth_factor = 2.0
	reject_large_step = True
    cutback_factor = 1.0
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
[]

