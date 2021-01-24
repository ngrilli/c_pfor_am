# Sandfeld and Zaiser model with zero curvature
# only first two equations with edge GND density are solved

# Stefan Sandfeld and Michael Zaiser
# Pattern formation in a minimal model of continuum
# dislocation plasticity
# Modelling Simul. Mater. Sci. Eng. 23 (2015) 065005 (18pp)

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 20
  ny = 20
  nz = 1
  xmax = 2.0
  ymax = 2.0
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
  
  [./rho_gnd_edge] # rho_x in the paper
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = FunctionIC
      function = init_rho_gnd_edge
    [../]
  [../]
  
  [./rho_t] # rho_t in the paper
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = FunctionIC
      function = init_rho_t
    [../]
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
  [./init_rho_gnd_edge]
    type = ParsedFunction
	value = 'if(0.999*x,0.0,1.0)+if(0.999*(2.0-x),0.0,-1.0)+if(x-1.0,0.0,2.0*(1.0-x))'
  [../]
  [./init_rho_t]
    type = ParsedFunction
	value = 'if(0.999*x,0.0,1.0)+if(0.999*(2.0-x),0.0,1.0)+if(x-1.0,0.0,2.0*abs(x-1.0))'
  [../]
  [./disp_load]
    type = ParsedFunction
    value = 'max(0.05*t,0.05*0.002)'
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
  [./drho_gnd_edge_dt]
    type = TimeDerivative
    variable = rho_gnd_edge
  [../]
  [./rho_gnd_edge_advection]
    type = ConservativeAdvectionCoupled
    variable = rho_gnd_edge
	rho_coupled = rho_t
    upwinding_type = none
	slip_sys_index = 0
	dislo_sign = positive
	dislo_character = edge
  [../]
  [./drho_t_dt]
    type = TimeDerivative
    variable = rho_t
  [../]
  [./rho_t_advection]
    type = ConservativeAdvectionCoupled
    variable = rho_t
	rho_coupled = rho_gnd_edge
    upwinding_type = none
	slip_sys_index = 0
	dislo_sign = positive
	dislo_character = edge
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
  
  [./rho_t_lr]
    type = DirichletBC
    variable = rho_t
    boundary = 'left right'
    value = 1.0  
  [../]
  
  [./rho_gnd_edge_left]
    type = DirichletBC
    variable = rho_gnd_edge
    boundary = left
    value = 1.0  
  [../]
  
  [./rho_gnd_edge_right]
    type = DirichletBC
    variable = rho_gnd_edge
    boundary = right
    value = -1.0  
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
    reference_temperature = '303.0'
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
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg          31'
  line_search = 'default'
  l_max_its = 50
  nl_max_its = 50
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  #l_tol = 1e-8

  start_time = 0.0
  end_time = 0.0002 #0.01
  dt = 0.0001
  dtmin = 0.00000001
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 2
  [../]
[]
