# Sandfeld and Zaiser model with constant curvature
# all equations solved

# Stefan Sandfeld and Michael Zaiser
# Pattern formation in a minimal model of continuum
# dislocation plasticity
# Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp)

# A2 matrix is taken from
# Katrin Schulz, Lydia Wagner, Christian Wieners
# A mesoscale continuum approach of dislocation dynamics and the
# approximation by a Runge-Kutta discontinuous Galerkin method
# International Journal of Plasticity 120 (2019) 248–261

# test with single dislocation loop

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 20
  ny = 20
  nz = 1
  xmax = 20.0
  ymax = 20.0
  zmax = 1.0
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
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = DislocationLoopsIC
      centrex = '10.0'
      centrey = '10.0'
      radii = '3.0'
      width = '1.0'
      rho_max = '0.01'
      variable_type = 'qtot'
	  slip_sys_index = 0
	  read_prop_user_object = prop_read
	  slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
	  nss = 12
    [../]
  [../]
  
  [./rho_t] # rho_t in the paper
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = DislocationLoopsIC
      centrex = '10.0'
      centrey = '10.0'
      radii = '3.0'
      width = '1.0'
      rho_max = '0.01'
      variable_type = 'rhotot'
	  slip_sys_index = 0
	  read_prop_user_object = prop_read
	  slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
	  nss = 12
    [../]
  [../]
  
  [./rho_gnd_edge] # rho_x in the paper
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = DislocationLoopsIC
      centrex = '10.0'
      centrey = '10.0'
      radii = '3.0'
      width = '1.0'
      rho_max = '0.01'
      variable_type = 'rhoedgegnd'
	  slip_sys_index = 0
	  read_prop_user_object = prop_read
	  slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
	  nss = 12
    [../]
  [../]
  
  [./rho_gnd_screw] # rho_x in the paper
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = DislocationLoopsIC
      centrex = '10.0'
      centrey = '10.0'
      radii = '3.0'
      width = '1.0'
      rho_max = '0.01'
      variable_type = 'rhoscrewgnd'
	  slip_sys_index = 0
	  read_prop_user_object = prop_read
	  slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
	  nss = 12
    [../]
  [../]
[]

[AuxVariables]
  [./temp]
    order = CONSTANT
    family = MONOMIAL
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

  [./gss1]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./temperature_load]
    type = ParsedFunction
    value = '303.0'
  [../]
  [./disp_load]
    type = ParsedFunction
    value = 'min(0.00005*t,0.0005)'
  [../]
  [./dts]
    type = PiecewiseConstant
    x = '0.0  100.0'
    y = '0.001  0.001'
  [../]
  [./rho_t_init]
    type = ParsedFunction
    value = '1.0+0.05*(cos(25*x)+cos(25*y))'
  [../]
  [./q_t_init]
    type = ParsedFunction
    value = '0.001+0.00005*(cos(25*x)+cos(25*y))'
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
  
  [./drho_t_dt]
    type = TimeDerivative
    variable = rho_t
  [../]
  [./rho_gnd_edge_advection]
    type = ConservativeAdvectionCoupled
    variable = rho_t
	rho_coupled = rho_gnd_edge
    upwinding_type = none
	slip_sys_index = 0
	dislo_sign = positive
	dislo_character = edge
  [../]
  [./rho_gnd_screw_advection]
    type = ConservativeAdvectionCoupled
    variable = rho_t
	rho_coupled = rho_gnd_screw
    upwinding_type = none
	slip_sys_index = 0
	dislo_sign = positive
	dislo_character = screw
  [../]  
  [./rho_t_curvature_mult]
    type = CurvatureMultiplication
    variable = rho_t
	curvature = q_t
	slip_sys_index = 0
  [../]
  [./rho_t_diffusion]
    type = CoefDiffusion
    variable = rho_t
    coef = 0.5
  [../]
  
  [./drho_gnd_edge_dt]
    type = TimeDerivative
	variable = rho_gnd_edge
  [../]
  [./rho_t_advection_edge]
    type = ConservativeAdvectionCoupled
    variable = rho_gnd_edge
	rho_coupled = rho_t
    upwinding_type = none
	slip_sys_index = 0
	dislo_sign = positive
	dislo_character = edge
  [../]
  [./rho_gnd_edge_diffusion]
    type = CoefDiffusion
    variable = rho_gnd_edge
    coef = 0.5
  [../]
  
  [./drho_gnd_screw_dt]
    type = TimeDerivative
	variable = rho_gnd_screw
  [../]
  [./rho_t_advection_screw]
    type = ConservativeAdvectionCoupled
    variable = rho_gnd_screw
	rho_coupled = rho_t
    upwinding_type = none
	slip_sys_index = 0
	dislo_sign = positive
	dislo_character = screw
  [../]
  [./rho_gnd_screw_diffusion]
    type = CoefDiffusion
    variable = rho_gnd_screw
    coef = 0.5
  [../]
  
  [./dq_t_dt]
    type = TimeDerivative
    variable = q_t
  [../]
  [./curvature_advection_edge]
    type = CurvatureAdvection
	variable = q_t
	rho_gnd = rho_gnd_edge
	rho_tot = rho_t
	slip_sys_index = 0
	dislo_character = edge
  [../]
  [./curvature_advection_screw]
    type = CurvatureAdvection
	variable = q_t
	rho_gnd = rho_gnd_screw
	rho_tot = rho_t
	slip_sys_index = 0
	dislo_character = screw
  [../]
  [./A2TraceX]
    type = A2Trace
	variable = q_t
    rho_tot = rho_t
	dv_dx = ddislovdx
	dv_dy = ddislovdy
	slip_sys_index = 0
	dislo_character = edge
  [../]
  [./A2TraceY]
    type = A2Trace
	variable = q_t
    rho_tot = rho_t
	dv_dx = ddislovdx
	dv_dy = ddislovdy
	slip_sys_index = 0
	dislo_character = screw
  [../]
  [./A2deviatoricX]
    type = A2Deviatoric
	variable = q_t
	rho_gnd_edge = rho_gnd_edge
	rho_gnd_screw = rho_gnd_screw
	dv_dx = ddislovdx
	dv_dy = ddislovdy
	slip_sys_index = 0
	dislo_character = edge	
  [../]
  [./A2deviatoricY]
    type = A2Deviatoric
	variable = q_t	
	rho_gnd_edge = rho_gnd_edge
	rho_gnd_screw = rho_gnd_screw
	dv_dx = ddislovdx
	dv_dy = ddislovdy
	slip_sys_index = 0	
	dislo_character = screw	
  [../]
  [./q_t_diffusion]
    type = CoefDiffusion
    variable = q_t
    coef = 0.5
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
    index = 0
    execute_on = timestep_end
  [../]
  
  [./ddislovdx]
    type = VariableGradientComponent
	variable = ddislovdx
	gradient_variable = dislov
	component = 'x'
  [../]
  
  [./ddislovdy]
    type = VariableGradientComponent
	variable = ddislovdy
	gradient_variable = dislov
	component = 'y'
  [../]

  [./stress_xz]
    type = RankTwoAux
    variable = stress_xz
    rank_two_tensor = stress
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  [../]

  [./gss1]
    type = MaterialStdVectorAux
    variable = gss1
    property = gss
    index = 0
    execute_on = timestep_end
  [../]

  [./fp_xz]
    type = RankTwoAux
    variable = fp_xz
    rank_two_tensor = fp
    index_j = 2
    index_i = 0
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
    type = FunctionDirichletBC
    variable = disp_x
    boundary = front
    function = disp_load
  [../]

  [./Periodic]

    # Can use auto_direction with Generated Meshes
    [./auto_rho_t_boundary]
      variable = rho_t
      auto_direction = 'x y'
    [../]
    [./auto_q_t_boundary]
      variable = q_t
      auto_direction = 'x y'
    [../]
    [./auto_rho_gnd_edge_boundary]
      variable = rho_gnd_edge
      auto_direction = 'x y'
    [../]
    [./auto_rho_gnd_screw_boundary]
      variable = rho_gnd_screw
      auto_direction = 'x y'
    [../]

  [../]
[]
 
[Postprocessors]

[]

[Materials]
  [./crysp]
    type = FiniteStrainCrystalPlasticityDislo
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 12 # Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' # slip rate equations parameters
	hprops = '1.0 3629.0 216.0 300.5 2.5' # hardening properties
    gprops = '1 12 216.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    thermal_expansion = '44.73e-6'
    reference_temperature = '303.0'
    temp = temp
    slip_incr_tol = 0.001
    maxiter = 200
    maxitergss = 200
    maximum_substep_iteration = 5
    gen_random_stress_flag = true
    #rtol = '1.0e-4'
    #abs_tol = '1.0e-4'
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
	burgers_vector_mag = 0.00025 # um
	shear_modulus_hardening = 10000.0 # MPa
    dislo_max_velocity = 1.0
	rho_edge_pos_1 = rho_t
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

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor -pc_hypre_boomeramg_print_statistics'

  petsc_options_value = 'hypre boomeramg 51 0.7 4 5 25 PMIS ext+i 2 0.3 0'

  line_search = 'default'

  automatic_scaling = true
  
  l_max_its = 50
  nl_max_its = 10
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  
  [./TimeStepper]
    type = FunctionDTGrowth
    function = dts
	cutback_factor_at_failure = 0.1
	growth_factor = 1.2
  [../]

  start_time = 0.0
  end_time = 0.001 #100.0
  
  dtmin = 1.0e-50
  timestep_tolerance = 1.0e-50
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 1 #30
  [../]
[]
