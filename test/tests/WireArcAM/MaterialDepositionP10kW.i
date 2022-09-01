# adding material on top of the substrate
# and use the ellipsoid heat source to simulate WAAM

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0.0
    xmax = 0.2
    ymin = 0.0
    ymax = 0.2
    zmin = 0.0
    zmax = 0.1
    nx = 80
    ny = 80
    nz = 40
    elem_type = HEX8
  [../]
  [./back_domain]
    input = gen
    type = SubdomainBoundingBoxGenerator
    bottom_left = '0.0 0.0 0.0'
    top_right = '0.2 0.2 0.05'
    block_id = 1
  [../]
  [./front_domain]
    input = back_domain
    type = SubdomainBoundingBoxGenerator
    bottom_left = '0.0 0.0 0.05'
    top_right = '0.2 0.2 0.1'
    block_id = 2
  [../]
  [./sidesets]
    input = front_domain
    type = SideSetsAroundSubdomainGenerator
    normal = '0 0 1'
    block = 1
    new_boundary = 'moving_interface'
  [../]  
[]

# the new elements will be initialized at temperature 300 K
[Variables]
  [./temp]
    initial_condition = 300.0
	block = '1'
  [../]
[]

# thermal conduction problem is solved only in the active domain
[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
	block = '1'
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    thermal_conductivity = thermal_conductivity
	block = '1'
  [../]
  [./heat_source]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
	block = '1'
  [../]
[]

[AuxVariables]
  [./level_set_var]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

# level set is initialized as 1.0
# only in the substrate
[ICs]
  [./temperature_ic_back]
    type = ConstantIC
    variable = level_set_var
    value = 1.0
	block = '1'
  [../]
  [./temperature_ic_front]
    type = ConstantIC
    variable = level_set_var
	value = 0.0
    block = '2'
  [../]
[]

[AuxKernels]
  [./add_material_in_the_ellipsoid]
    type = FunctionPathEllipsoidAux
	variable = level_set_var
	level_set_var = level_set_var
	rx = 0.01
    ry = 0.01
    rz = 0.01
    function_x= path_x
    function_y= path_y
    function_z= path_z
	level_set_activation_threshold = 0.3
  [../]
[]

# Neumann BC
# no heat exiting the volume
[BCs]
[]

# apply the ellipsoid heat source
# on the same path that is used
# for adding elements
[Materials]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 420 # J/KgK steel
    thermal_conductivity = 45.0 # W/mK steel
  [../]
  [./density]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = 7850.0 # Kg/m^3 steel
  [../]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
	rx = 0.01
    ry = 0.01
    rz = 0.01
    power = 10000 # W -> CO2 laser
    efficiency = 1.0
    factor = 1
    function_x= path_x
    function_y= path_y
    function_z= path_z
  [../]
[]

[Functions]
  [./path_x]
    type = PiecewiseLinear
    x = '0.00 10.0 11.0 21.0 22.0 32.0 33.0 35.0 45.0 46.0 56.0 57.0 67.0'
    y = '0.00 0.20 0.20 0.00 0.00 0.20 0.20 0.20 0.00 0.00 0.20 0.20 0.00'
  [../]
  [./path_y]
    type = PiecewiseLinear
    x = '0.00 10.0 11.0 21.0 22.0 32.0 33.0 35.0 45.0 46.0 56.0 57.0 67.0'
    y = '0.075 0.075 0.10 0.10 0.125 0.125 0.125 0.075 0.075 0.10 0.10 0.125 0.125'
  [../]
  [./path_z]
    type = PiecewiseLinear
    x = '0.00 10.0 11.0 21.0 22.0 32.0 33.0 35.0 45.0 46.0 56.0 57.0 67.0'
    y = '0.05 0.05 0.05 0.05 0.05 0.05 0.075 0.075 0.075 0.075 0.075 0.075 0.075'
  [../]
[]

# User Object to activate elements
[UserObjects]
  [./activated_elem_uo]
    type = ActDeactElementsCoupled
    execute_on = timestep_begin
    coupled_var = level_set_var
    activate_value = 0.5
	activate_type = 'above'
    active_subdomain_id = 1
	deactive_subdomain_id = 2
    expand_boundary_name = 'moving_interface'
  [../]
[]

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor -pc_hypre_boomeramg_print_statistics'
  petsc_options_value = 'hypre boomeramg 51 0.7 4 5 25 PMIS ext+i 2 0.3 0'

  l_max_its = 50
  nl_max_its = 50

  end_time = 67
  dt = 0.1
  dtmin = 1e-20
[]

[Outputs]
  print_linear_residuals = true
  perf_graph = true
  [./exodus]
    type = Exodus
    interval = 10
  [../]
[]
