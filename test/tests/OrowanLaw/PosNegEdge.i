[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 20
  ny = 10
  nz = 1
  xmax = 2.0
  ymax = 1.0
  zmax = 0.1
  elem_type = HEX8
  displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
  
  # Dislocation densities
  
  [./rho_edge_pos_12] 
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = FunctionIC
      function = init_rho_edge_pos
    [../]
  [../]
  
  [./rho_edge_neg_12] 
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = FunctionIC
      function = init_rho_edge_neg
    [../]
  [../]
[]

[AuxVariables]
  [./temp]
    order = FIRST
    family = LAGRANGE
  [../]

  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./gss]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./dislov]
    order = CONSTANT
    family = MONOMIAL
  [../]  
  
  [./slip_incr_out]
    order = CONSTANT
    family = MONOMIAL
  [../]  
[]

[Functions]
  [./temperature_load]
    type = ParsedFunction
    value = '298.0'
  [../]
  [./init_rho_edge_pos]
    type = ParsedFunction
	value = 'if(3.0*(x-1.05),0.0,1.0)*if(3.0*(y-0.5),0.0,1.0)'
  [../]
  [./init_rho_edge_neg]
    type = ParsedFunction
	value = 'if(3.0*(x-0.95),0.0,1.0)*if(3.0*(y-0.5),0.0,1.0)'
  [../]
  [./disp_load]
    type = ParsedFunction
    value = '0.1*t'
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
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    add_variables = true
  [../]
  [./drho_pos_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_12
  [../]
  [./rho_pos_advection]
    type = ConservativeAdvectionSchmid
    variable = rho_edge_pos_12
    upwinding_type = full
	slip_sys_index = 11
	dislo_sign = positive
  [../]
  [./drho_neg_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_12
  [../]
  [./rho_neg_advection]
    type = ConservativeAdvectionSchmid
    variable = rho_edge_neg_12
    upwinding_type = full
	slip_sys_index = 11
	dislo_sign = negative
  [../]
[]

[AuxKernels]

  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
    index_j = 1
    index_i = 0
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]  

  [./fp_xy]
    type = RankTwoAux
    variable = fp_xy
    rank_two_tensor = fp
    index_j = 1
    index_i = 0
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./tempfuncaux]
    type = FunctionAux
    variable = temp
    function = temperature_load
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./gss]
    type = MaterialStdVectorAux
    variable = gss
    property = gss
    index = 11
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]
  
  [./dislov]
    type = MaterialStdVectorAux
    variable = dislov
    property = dislo_velocity
    index = 11
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]
  
  [./slip_incr_out]
    type = MaterialStdVectorAux
    variable = slip_incr_out
    property = slip_incr_out
    index = 11
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

[]

[BCs]

  [./z_back]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./y_bott]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  
  [./x_bott]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0.0
  [../]
  
  [./disp_top]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = top
    function = disp_load
  [../]

[]
 
[Postprocessors]

[]

[Materials]
  [./crysp]
    type = FiniteStrainCrystalPlasticityDislo
    block = 0
    gtol = 1e-2
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 12 # Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' # slip rate equations parameters
	hprops = '1.0 3629.0 216.0 300.5 2.5' # hardening properties
    gprops = '1 12 216.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    thermal_expansion = '17e-6'
    reference_temperature = '298.0'
    temp = temp
# Calibrated using table 1 in:
# M.R. DAYMOND and P.J. BOUCHARD
# Elastoplastic Deformation of 316 Stainless Steel Under
# Tensile Loading at Elevated Temperatures
# METALLURGICAL AND MATERIALS TRANSACTIONS A
# VOLUME 37A, JUNE 2006—1873
	dCRSS_dT_A = 0.53
	dCRSS_dT_B = 0.47
	dCRSS_dT_C = 0.008
# Hull, Bacon, Dislocations, figure 3.11a
	dislo_mobility = 1.0 # um/s/MPa
	burgers_vector_mag = 0.001
	rho_edge_pos_12 = rho_edge_pos_12
	rho_edge_neg_12 = rho_edge_neg_12
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
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
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
  end_time = 0.002 # 0.35
  dt = 0.001
  dtmin = 0.00001
[]

[Outputs]
  [./out]
    type = Exodus
    time_step_interval = 2
  [../]
[]
