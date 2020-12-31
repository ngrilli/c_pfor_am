[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim =  3
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 1.0
    nx= 2
    ny= 1
    nz= 1
    elem_type = HEX8
  [../]
  [./active_domain]
    input = gen
    type = SubdomainBoundingBoxGenerator
    bottom_left = '0.0 0.0 0.0'
    top_right = '1.0 1.0 1.0'
    block_id = 1
  [../]
  [./deactivated_domain]
    input = active_domain
    type = SubdomainBoundingBoxGenerator
    bottom_left = '1.0 0.0 0.0'
    top_right = '2.0 1.0 1.0'
    block_id = 2
  [../]
  [./sidesets]
    input = deactivated_domain
    type = SideSetsAroundSubdomainGenerator
    normal = '1 0 0'
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
  [./temp]
    order = CONSTANT
    family = MONOMIAL
	[./InitialCondition]
      type = ConstantIC
      value = 1000.0
    [../]
  [../]

  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
[]

[Functions]
  [./disp_load]
    type = PiecewiseLinear
    x = '0.0 1.0 8.0'
    y = '0.0 0.0001 0.0001'
  [../]
[]

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'euler_ang_test.inp'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    ngrain = 1
    read_type = indexgrain
  [../]
  [./temperature_read]
    type = LaserTempReadFile
	temperature_file_name = 'temperature1El.txt'
	temperature_num_step = 9
  [../]
  [./activated_elem_uo]
    type = ActDeactElementsMelting
    execute_on = timestep_begin
    temperature = temp
    active_subdomain_id = 1
	deactive_subdomain_id = 2
    expand_boundary_name = 'moving_interface'
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

  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]

  [./tempfuncaux]
    type = TempActDeactElemsAux
    variable = temp
    temperature_read_user_object = temperature_read
	temperature_time_step = 1.0
  [../]
[]

[BCs]

  [./z_back]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./y_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./x_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  
  [./x_right]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0.0
  [../]
  
  [./z_load]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = disp_load
  [../]

[]
 
[Postprocessors]

[]

[Materials]
  [./crysp]
    type = FiniteStrainCrystalPlasticityThermal
    block = '1'
    gtol = 1e-2
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 12 # Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' # slip rate equations parameters 
# Calibrated by comparing with Fig 2b in:
# Wen Chen et al.
# Microscale residual stresses in additively
# manufactured stainless steel
# NATURE COMMUNICATIONS (2019) 10:4338
    hprops = '1.0 3839.0 213.0 302.0 2.5' # hardening properties
    gprops = '1 12 213.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    thermal_expansion = '0.0e-6'
    reference_temperature = '298.0'
    temp = temp
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
    dCTE_dT='0.0'
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
    temp = temp
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

[Preconditioning]
  active = 'smp'
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
  petsc_options_value = 'hypre    boomeramg          31'
  line_search = 'none'
  l_max_its = 50
  nl_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  l_tol = 1e-8

  start_time = 0.0
  end_time = 0.1 #8.0
  dt = 0.1
  dtmin = 0.01
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 1 #2
  [../]
[]
