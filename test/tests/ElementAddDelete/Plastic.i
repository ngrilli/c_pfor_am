[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim =  3
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 1.0
    nx= 10
    ny= 10
    nz= 10
  [../]
  [./left_domain]
    input = gen
    type = SubdomainBoundingBoxGenerator
    bottom_left = '0.0 0.0 0.0'
    top_right = '1.0 1.0 0.5'
    block_id = 1
  [../]
  [./right_domain]
    input = left_domain
    type = SubdomainBoundingBoxGenerator
    bottom_left = '0.0 0.0 0.5'
    top_right = '1.0 1.0 1.0'
    block_id = 2
  [../]
  [./sidesets]
    input = right_domain
    type = SideSetsAroundSubdomainGenerator
    normal = '0 0 1'
    block = 1
    new_boundary = 'moving_interface'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    # Use the initial Condition block underneath the variable
    # for which we want to apply this initial condition
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
    block = '1'
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]

  [./disp_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
[]

[AuxVariables]
  [./temperature]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
	block = '1'
  [../]
  
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
	block = '1'
  [../]
  
  [./fp_xx]
    order = CONSTANT
    family = MONOMIAL
	block = '1'
  [../]
[]

[ICs]
  [./temperature_ic]
    type = FunctionIC
    variable = temperature
    function = temperature_load_init
  [../]
[]

[Functions]
  [./temperature_load]
    type = ParsedFunction
    value = 'if((t-100.0)/200.0,0.0,303.0+100.0*z)+if((t-250.0)/100.1,0.0,303.0+(t-100.0)*z)+if((t-350.0)/99.9,0.0,303.0+(500.0-t)*z)+if((t-450.0)/100.1,0.0,303.0+100.0*z)'
  [../]
  [./temperature_load_init]
    type = ParsedFunction
    value = '303.0+(100.0*z)'
  [../]
  [./disp_load]
    type = PiecewiseLinear
    x = '0.0 100.0 200.0 300.0 400.0 500.0'
    y = '0.0 0.01 0.0 0.0 0.0 0.01'
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    add_variables = true
	block = '1'
  [../]
[]

[AuxKernels]
  [./tempfuncaux]
    type = FunctionAux
    variable = temperature
    function = temperature_load
  [../]
  
  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
	block = '1'
  [../]
  
  [./strain_xx]
    type = RankTwoAux
    variable = strain_xx
    rank_two_tensor = total_strain
    index_j = 0
    index_i = 0
    execute_on = timestep_end
	block = '1'
  [../]
  
  [./fp_xx]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = fp
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./crysp]
    type = FiniteStrainCrystalPlasticityThermal
    gtol = 1e-2
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 12 #Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' # slip rate equations parameters
    hprops = '1.0 3839.0 213.0 302.0 2.5' # hardening properties
    gprops = '1 12 213.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    #thermal_expansion = '44.73e-6'
    #reference_temperature = '293.0'
    #temp = temperature
    slip_incr_tol = 0.001
    maxiter = 250
    maxitergss = 250
    maximum_substep_iteration = 6
    gen_random_stress_flag = true
# Calibrated using table 1 in:
# M.R. DAYMOND and P.J. BOUCHARD
# Elastoplastic Deformation of 316 Stainless Steel Under
# Tensile Loading at Elevated Temperatures
# METALLURGICAL AND MATERIALS TRANSACTIONS A
# VOLUME 37A, JUNE 2006—1873
    dCRSS_dT_A = 0.53
    dCRSS_dT_B = 0.47
    dCRSS_dT_C = 0.008
# Calibrated using table 1 in:
# W.Jiang, Y.Zhang and W.Woo
# Using heat sink technology to decrease residual stress
# in 316L stainless steel welding joint:Finite element simulation
# Int.J.Press.Vessel.Pip.
# VOLUME 92, pp.56-62, 2012
    #dCTE_dT='0.01011e-6'
	block = '1'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPGrain
# Elastic constants of 316L SS from:
# Clausen, B., Lorentzen, T. and Leffers, T.
# Self-consistent modelling of the plastic
# deformation of FCC polycrystals and its implications for diffraction
# measurements of internal stresses.
# Acta Mater. 46, 3087–3098 (1998).
    C_ijkl = '2.046e5 1.377e5 1.377e5 2.046e5 1.377e5 2.046e5 1.262e5 1.262e5 1.262e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
    #temp = temp
    dC11_dT = 0.0004415
    dC12_dT = 0.0003275
    dC44_dT = 0.0004103
	block = '1'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
	block = '1'
  [../]

  [./dummy_mat_inactive]
    type = GenericConstantMaterial
    prop_names = 'dummy_mat'
    prop_values = '0.0'
	block = '2'
  [../]
[]

[BCs]
  [./z0]
    type = DirichletBC
    variable = disp_z
    boundary = back
	value = 0.0
  [../]

  [./y0]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./x0]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  
  [./x1]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = disp_load
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
  solve_type = 'PJFNK'
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg          101'
  line_search = 'none'
  l_max_its = 50
  nl_max_its = 20
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  l_tol = 1e-8

  start_time = 0.0
  end_time = 1.0 #500.0
  dt = 1.0
  dtmin = 1e-4
[]

[UserObjects]
  [./activated_elem_uo]
    type = ActDeactElementsCoupled
    execute_on = timestep_begin
    coupled_var = temperature
    activate_value = 403.0
	activate_type = below
    active_subdomain_id = 1
	deactive_subdomain_id = 2
    expand_boundary_name = 'moving_interface'
  [../]
  
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'euler_ang_test.inp'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    ngrain = 1
    read_type = indexgrain
  [../]
[]

[Outputs]
  exodus = true
  interval = 1 #10
[]

