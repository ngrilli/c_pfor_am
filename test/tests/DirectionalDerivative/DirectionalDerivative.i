# Calculate directional derivative of a given function

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmax = 1.0
    ymax = 1.0
    zmax = 1.0
    elem_type = HEX8
	displacements = 'disp_x disp_y disp_z'
  [../]
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
  
  [./var_for_grad] # variable to take the gradient
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = 'x+y+z'
    [../]
  [../]
[]

[AuxVariables]
  [./temp]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./dvar_dedge]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  
  [./dvar_dscrew]
    order = CONSTANT
    family = MONOMIAL  
  [../]
[]

[Functions]
  [./temperature_load]
    type = ParsedFunction
    value = '303.0'
  [../]
  [./init_var_for_grad]
    type = ParsedFunction
    value = 'min(0.0001*t,0.0002)'
  [../]
  [./dts]
    type = PiecewiseConstant
    x = '0.0 1.0'
    y = '0.5 0.5'
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
    variable = var_for_grad
  [../]
[]

[AuxKernels]
  [./tempfuncaux]
    type = FunctionAux
    variable = temp
    function = temperature_load
  [../]
  
  [./dvar_dedge]
    type = DirectionalDerivative
	variable = dvar_dedge
	gradient_variable = var_for_grad
	slip_sys_index = 0
	dislo_character = edge
  [../]
  
  [./dvar_dscrew]
    type = DirectionalDerivative
	variable = dvar_dscrew
    gradient_variable = var_for_grad
	slip_sys_index = 0
	dislo_character = screw
  [../]
[]

[BCs]
  [./x_bc]
    type = DirichletBC
    variable = disp_x
    boundary = 'top bottom left right front back'
    value = 0.0 
  [../]
  [./y_bc]
    type = DirichletBC
    variable = disp_y
    boundary = 'top bottom left right front back'
    value = 0.0 
  [../]
  [./z_bc]
    type = DirichletBC
    variable = disp_z
    boundary = 'top bottom left right front back'
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
	dislo_mobility = 1 # um/s/MPa
	burgers_vector_mag = 0.00025 # um
	shear_modulus_hardening = 10000.0 # MPa
    dislo_max_velocity = 1.0
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
  nl_max_its = 50
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8
  
  [./TimeStepper]
    type = FunctionDTGrowth
    function = dts
	cutback_factor_at_failure = 0.1
	growth_factor = 1.2
  [../]

  start_time = 0.0
  end_time = 1.0
  
  dtmin = 1.0e-50
  timestep_tolerance = 1.0e-50
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 1
  [../]
[]
