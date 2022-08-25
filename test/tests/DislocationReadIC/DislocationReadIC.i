# four variables dislocation transport model
# read initial values of the variables from file

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 5
  nz = 1
  xmax = 0.5
  ymax = 0.5
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
  
  [./q_t] # curvature density q_t in the paper
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = DislocationReadIC
      read_dislocation_user_object = dislocation_read
      variable_type = 'qtot'
    [../]
  [../]
  
  [./rho_t] # rho_t in the paper
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = DislocationReadIC
      read_dislocation_user_object = dislocation_read
      variable_type = 'rhotot'
    [../]
  [../]
  
  [./rho_gnd_edge] # rho_x in the paper
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = DislocationReadIC
      read_dislocation_user_object = dislocation_read
      variable_type = 'rhoedgegnd'
    [../]
  [../]
  
  [./rho_gnd_screw] # rho_x in the paper
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = DislocationReadIC
      read_dislocation_user_object = dislocation_read
      variable_type = 'rhoscrewgnd'
    [../]
  [../]
[]

[AuxVariables]
  [./temp]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = ConstantIC
      value = '303.0'
    [../]
  [../]
  
  [./dislov]
    order = FIRST
    family = MONOMIAL
  [../] 
  [./ddislovdx]
    order = CONSTANT
    family = MONOMIAL
  [../] 
  [./ddislovdy]
    order = CONSTANT
    family = MONOMIAL
  [../] 
  
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./gss8]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./temperature_load]
    type = ParsedFunction
    value = '303.0'
  [../]
  [./dts]
    type = PiecewiseConstant
    x = '0.0 20.0'
    y = '1.0 1.0'
  [../]
  [./velocity_x]
    type = ParsedFunction
    value = '0.0000004'
  [../]
[]

[UserObjects]
  [./dislocation_read]
    type = PropertyReadFile
    prop_file_name = 'initial_dislocations.txt'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 4
    read_type = element
  [../]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'euler_angles.txt'
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
  
  [./drho_t_dt]
    type = TimeDerivative
    variable = rho_t
  [../]
  [./rho_gnd_edge_advection]
    type = ConservativeAdvectionCoupled
    variable = rho_t
	rho_coupled = rho_gnd_edge
	use_displaced_mesh = false
    upwinding_type = none
	slip_sys_index = 7
	dislo_sign = positive
	dislo_character = edge
  [../]
  [./rho_gnd_screw_advection]
    type = ConservativeAdvectionCoupled
    variable = rho_t
	rho_coupled = rho_gnd_screw
	use_displaced_mesh = false
    upwinding_type = none
	slip_sys_index = 7
	dislo_sign = positive
	dislo_character = screw
  [../]  

  [./drho_gnd_edge_dt]
    type = TimeDerivative
	variable = rho_gnd_edge
  [../]
  [./rho_t_advection_edge]
    type = ConservativeAdvectionCoupled
    variable = rho_gnd_edge
	rho_coupled = rho_t
	use_displaced_mesh = false
    upwinding_type = none
	slip_sys_index = 7
	dislo_sign = positive
	dislo_character = edge
  [../]

  [./drho_gnd_screw_dt]
    type = TimeDerivative
	variable = rho_gnd_screw
  [../]
  [./rho_t_advection_screw]
    type = ConservativeAdvectionCoupled
    variable = rho_gnd_screw
	rho_coupled = rho_t
	use_displaced_mesh = false
    upwinding_type = none
	slip_sys_index = 7
	dislo_sign = positive
	dislo_character = screw
  [../]

  [./dq_t_dt]
    type = TimeDerivative
    variable = q_t
	use_displaced_mesh = false
  [../]
  [./curvature_advection_edge]
    type = CurvatureAdvection
	variable = q_t
	use_displaced_mesh = false
	rho_gnd = rho_gnd_edge
	rho_tot = rho_t
	slip_sys_index = 7
	dislo_character = edge
    rho_tot_tol = 0.2
  [../]
  [./curvature_advection_screw]
    type = CurvatureAdvection
	variable = q_t
	use_displaced_mesh = false
	rho_gnd = rho_gnd_screw
	rho_tot = rho_t
	slip_sys_index = 7
	dislo_character = screw
    rho_tot_tol = 0.2
  [../]

  [./q_t_diffusion]
    type = CoefDiffusion
    variable = q_t
    coef = 0.001
  [../]
  [./rho_t_diffusion]
    type = CoefDiffusion
    variable = rho_t
     use_displaced_mesh = false
    coef = 0.001
  [../]
  [./rho_gnd_edge_diffusion]
    type = CoefDiffusion
    variable = rho_gnd_edge
	use_displaced_mesh = false
    coef = 0.001
  [../]
  [./rho_gnd_screw_diffusion]
    type = CoefDiffusion
    variable = rho_gnd_screw
	use_displaced_mesh = false
    coef = 0.001
  [../]

[]

[AuxKernels]
  [./tempfuncaux]
    type = FunctionAux
    variable = temp
    function = temperature_load
  [../]
  
  [./dislov]
    type = MaterialStdVectorAux
    variable = dislov
    property = dislo_velocity
    index = 7
    execute_on = timestep_end
  [../]
  
  [./ddislovdx]
    type = DirectionalDerivative
    variable = ddislovdx
    gradient_variable = dislov
    slip_sys_index = 7
    dislo_character = edge
  [../]
 
  [./ddislovdy]
    type = DirectionalDerivative
    variable = ddislovdy
    gradient_variable = dislov
    slip_sys_index = 7
    dislo_character = screw
  [../]

  [./stress_xz]
    type = RankTwoAux
    variable = stress_xz
    rank_two_tensor = stress
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  [../]
  
  [./gss8]
    type = MaterialStdVectorAux
    variable = gss8
    property = gss
    index = 7
    execute_on = timestep_end
  [../]

[]

[BCs]
  [./z_back]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  [./x_back]
    type = DirichletBC
    variable = disp_x
    boundary = back
    value = 0.0
  [../]
  [./y_back]
    type = DirichletBC
    variable = disp_y
    boundary = back
    value = 0.0
  [../]
  
  [./z_front]
    type = DirichletBC
    variable = disp_z
    boundary = front
    value = 0.0
  [../]
  [./y_front]
    type = DirichletBC
    variable = disp_y
    boundary = front
    value = 0.0
  [../]
  [./x_front]
    type = PresetVelocity
    variable = disp_x
    boundary = front
    function = velocity_x
	use_displaced_mesh = false
  [../]
  
  [./Periodic]
    [./all_rho_t]
      variable = rho_t
      auto_direction = 'x y z'
    [../]
    [./all_rho_gnd_edge]
      variable = rho_gnd_edge
      auto_direction = 'x y z'
    [../]
    [./all_rho_gnd_screw]
      variable = rho_gnd_screw
      auto_direction = 'x y z'
    [../]
    [./all_q_t]
      variable = q_t
      auto_direction = 'x y z'
    [../]	
  [../] 
[]
 
[Postprocessors]

[]

[Materials]
  [./strain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./crysp]
    type = FiniteStrainCrystalPlasticityDislo
    slip_sys_file_name = 'input_slip_sys.txt' # no need to normalize vectors
    nss = 12 # Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 12 0.001 0.1' # slip rate equations parameters
	hprops = '1.0 3629.0 216.0 300.5 2.5' # hardening properties
    gprops = '1 12 216.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    thermal_expansion = '44.73e-6'
    reference_temperature = '303.0'
    temp = temp
    slip_incr_tol = 0.001
    maxiter = 200
    maxitergss = 200
    maximum_substep_iteration = 1
    gen_random_stress_flag = true
# Calibrated using table 1 in:
# M.R. DAYMOND and P.J. BOUCHARD
# Elastoplastic Deformation of 316 Stainless Steel Under
# Tensile Loading at Elevated Temperatures
# METALLURGICAL AND MATERIALS TRANSACTIONS A
# VOLUME 37A, JUNE 2006\E2\80?873
    dCRSS_dT_A = 0.53
    dCRSS_dT_B = 0.47
    dCRSS_dT_C = 0.008
# Hull, Bacon, Dislocations, figure 3.11a
    dislo_mobility = 0.000256 # um/us/MPa
    reduced_mobility = 0.00000001 # um/us/MPa
    burgers_vector_mag = 0.00025 # um
    shear_modulus_hardening = 86000.0 # MPa
    dislo_max_velocity = 0.003 # um/us
    rho_edge_pos_8 = rho_t
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPGrain
# Elastic constants of 316L SS from:
# Clausen, B., Lorentzen, T. and Leffers, T. 
# Self-consistent modelling of the plastic
# deformation of FCC polycrystals and its implications for diffraction
# measurements of internal stresses.
# Acta Mater. 46, 3087\E2\80?098 (1998).
    C_ijkl = '2.046e5 1.377e5 1.377e5 2.046e5 1.377e5 2.046e5 1.262e5 1.262e5 1.262e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
    temp = temp
    dC11_dT = 0.0004415
    dC12_dT = 0.0003275
    dC44_dT = 0.0004103
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
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor -pc_hypre_boomeramg_print_statistics'
  petsc_options_value = 'hypre boomeramg 51 0.7 4 5 25 PMIS ext+i 2 0.3 0'
  
  line_search = 'none'
  automatic_scaling = true
  
  l_max_its = 40
  nl_max_its = 40
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-5
  
  [./TimeStepper]
    type = FunctionDTGrowth
    function = dts
	cutback_factor_at_failure = 0.1
	growth_factor = 1.2
  [../]

  start_time = 0.0
  end_time = 1.0
  
  dtmin = 1.0e-20
  timestep_tolerance = 1.0e-20
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 1
  [../]
[]
