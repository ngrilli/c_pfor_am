# Sandfeld and Zaiser model with constant curvature
# all equations solved
# only curvature diffusion is missing

# Stefan Sandfeld and Michael Zaiser
# Pattern formation in a minimal model of continuum
# dislocation plasticity
# Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp)

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 40
  ny = 40
  nz = 1
  xmax = 2.0
  ymax = 2.0
  zmax = 0.05
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
      type = RandomIC
	  min = 0.05
	  max = 0.1
    [../]
  [../]
  
  [./rho_t] # rho_t in the paper
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = RandomIC
	  min = 0.5
	  max = 1.0
    [../]
  [../]
  
  [./rho_gnd_edge] # rho_x in the paper
    order = FIRST
    family = LAGRANGE
  [../]
  
  [./rho_gnd_screw] # rho_x in the paper
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./temp]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./dislov]
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
    value = 'min(0.000025*t,0.000025)'
  [../]
  [./dts]
    type = PiecewiseConstant
    x = '0.0  100.0'
    y = '0.0025  0.0025'
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
  
  [./rho_t_boundary]
    type = DirichletBC
    variable = rho_t
    boundary = 'bottom top left right'
	value = 0.75
  [../]
  
  [./q_t_boundary]
    type = DirichletBC
    variable = q_t
    boundary = 'bottom top left right'
	value = 0.0
  [../]
  
  [./rho_gnd_edge_boundary]
    type = DirichletBC
    variable = rho_gnd_edge
    boundary = 'bottom top left right'
	value = 0.0
  [../]

  [./rho_gnd_screw_boundary]
    type = DirichletBC
    variable = rho_gnd_screw
    boundary = 'bottom top left right'
	value = 0.0
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
    rtol = '1.0e-4'
    abs_tol = '1.0e-4'
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
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  
  line_search = 'default'
  automatic_scaling = true
  
  l_max_its = 50
  nl_max_its = 20
  nl_rel_tol = 5e-6
  nl_abs_tol = 1e-8
  
  [./TimeStepper]
    type = FunctionDTGrowth
    function = dts
	cutback_factor_at_failure = 0.1
	growth_factor = 1.2
  [../]

  start_time = 0.0
  end_time = 10.0
  
  dtmin = 0.00000001
  timestep_tolerance = 0.00000001
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 100
  [../]
[]
