[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmax = 1.0
  ymax = 1.0
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
  
  [./rho_edge_pos_1] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_2] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_3] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_4] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_5] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_6] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_7] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_8] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_9] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_10] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_11] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_pos_12] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_1] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_2] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_3] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_4] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_5] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_6] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_7] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_8] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_9] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_10] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_11] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_edge_neg_12] # Edge dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_1] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_2] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_3] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_4] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_5] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_6] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_7] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_8] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_9] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_10] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_11] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_pos_12] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_1] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_2] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_3] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_4] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_5] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_6] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_7] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_8] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_9] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_10] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_11] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  
  [./rho_screw_neg_12] # Screw dislocation density
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
[]

[AuxVariables]
  [./temp]
    order = FIRST
    family = LAGRANGE
  [../]
  
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./gss]
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
    value = '293.0'
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
  
  [./drho_edge_pos_1_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_1
  [../]
  
  [./drho_edge_pos_2_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_2
  [../]
  
  [./drho_edge_pos_3_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_3
  [../]
  
  [./drho_edge_pos_4_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_4
  [../]
  
  [./drho_edge_pos_5_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_5
  [../]
  
  [./drho_edge_pos_6_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_6
  [../]
  
  [./drho_edge_pos_7_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_7
  [../]
  
  [./drho_edge_pos_8_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_8
  [../]
  
  [./drho_edge_pos_9_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_9
  [../]
  
  [./drho_edge_pos_10_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_10
  [../]
  
  [./drho_edge_pos_11_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_11
  [../]
  
  [./drho_edge_pos_12_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_12
  [../]
  
  [./drho_edge_neg_1_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_1
  [../]
  
  [./drho_edge_neg_2_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_2
  [../]
  
  [./drho_edge_neg_3_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_3
  [../]
  
  [./drho_edge_neg_4_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_4
  [../]
  
  [./drho_edge_neg_5_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_5
  [../]
  
  [./drho_edge_neg_6_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_6
  [../]
  
  [./drho_edge_neg_7_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_7
  [../]
  
  [./drho_edge_neg_8_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_8
  [../]
  
  [./drho_edge_neg_9_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_9
  [../]
  
  [./drho_edge_neg_10_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_10
  [../]
  
  [./drho_edge_neg_11_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_11
  [../]
  
  [./drho_edge_neg_12_dt]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_12
  [../]

  [./drho_screw_pos_1_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_1
  [../]
  
  [./drho_screw_pos_2_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_2
  [../]
  
  [./drho_screw_pos_3_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_3
  [../]
  
  [./drho_screw_pos_4_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_4
  [../]
  
  [./drho_screw_pos_5_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_5
  [../]
  
  [./drho_screw_pos_6_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_6
  [../]
  
  [./drho_screw_pos_7_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_7
  [../]
  
  [./drho_screw_pos_8_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_8
  [../]
  
  [./drho_screw_pos_9_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_9
  [../]
  
  [./drho_screw_pos_10_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_10
  [../]
  
  [./drho_screw_pos_11_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_11
  [../]
  
  [./drho_screw_pos_12_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_pos_12
  [../]
  
  [./drho_screw_neg_1_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_1
  [../]
  
  [./drho_screw_neg_2_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_2
  [../]
  
  [./drho_screw_neg_3_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_3
  [../]
  
  [./drho_screw_neg_4_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_4
  [../]
  
  [./drho_screw_neg_5_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_5
  [../]
  
  [./drho_screw_neg_6_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_6
  [../]
  
  [./drho_screw_neg_7_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_7
  [../]
  
  [./drho_screw_neg_8_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_8
  [../]
  
  [./drho_screw_neg_9_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_9
  [../]
  
  [./drho_screw_neg_10_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_10
  [../]
  
  [./drho_screw_neg_11_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_11
  [../]
  
  [./drho_screw_neg_12_dt]
    type = MassLumpedTimeDerivative
    variable = rho_screw_neg_12
  [../]

[]

[AuxKernels]
  [./tempfuncaux]
    type = FunctionAux
    variable = temp
    function = temperature_load
    block = 'ANY_BLOCK_ID 0'
  [../]
  
  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  
  [./gss]
    type = MaterialStdVectorAux
    variable = gss
    property = gss
    index = 0
    execute_on = timestep_end
  [../]
  
  [./slip_incr_out]
    type = MaterialStdVectorAux
    variable = slip_incr_out
    property = slip_incr_out
    index = 0
    execute_on = timestep_end
  [../]
  
[]

[BCs]
  [./z0_back]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./y0_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./x0_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  [./z1_disp]
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
    type = FiniteStrainCrystalPlasticityDislo
    block = 0
    gtol = 2e-2
	slip_incr_tol = 5e-2
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 12 # Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' # slip rate equations parameters
	hprops = '1.0 3629.0 216.0 300.5 2.5' # hardening properties
    gprops = '1 12 216.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    thermal_expansion = '17e-6'
    reference_temperature = '293.0'
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
# Hull, Bacon, Dislocations, figure 3.11a
	dislo_mobility = 500.0 # um/s/MPa
	burgers_vector_mag = 0.000258 # um
	shear_modulus_hardening = 86000.0 # MPa
	rho_edge_pos_1 = rho_edge_pos_1
	rho_edge_pos_2 = rho_edge_pos_2
	rho_edge_pos_3 = rho_edge_pos_3
	rho_edge_pos_4 = rho_edge_pos_4
	rho_edge_pos_5 = rho_edge_pos_5
	rho_edge_pos_6 = rho_edge_pos_6
	rho_edge_pos_7 = rho_edge_pos_7
	rho_edge_pos_8 = rho_edge_pos_8
	rho_edge_pos_9 = rho_edge_pos_9
	rho_edge_pos_10 = rho_edge_pos_10
	rho_edge_pos_11 = rho_edge_pos_11
	rho_edge_pos_12 = rho_edge_pos_12
	rho_edge_neg_1 = rho_edge_neg_1
	rho_edge_neg_2 = rho_edge_neg_2
	rho_edge_neg_3 = rho_edge_neg_3
	rho_edge_neg_4 = rho_edge_neg_4
	rho_edge_neg_5 = rho_edge_neg_5
	rho_edge_neg_6 = rho_edge_neg_6
	rho_edge_neg_7 = rho_edge_neg_7
	rho_edge_neg_8 = rho_edge_neg_8
	rho_edge_neg_9 = rho_edge_neg_9
	rho_edge_neg_10 = rho_edge_neg_10
	rho_edge_neg_11 = rho_edge_neg_11
	rho_edge_neg_12 = rho_edge_neg_12
	rho_screw_pos_1 = rho_screw_pos_1
	rho_screw_pos_2 = rho_screw_pos_2
	rho_screw_pos_3 = rho_screw_pos_3
	rho_screw_pos_4 = rho_screw_pos_4
	rho_screw_pos_5 = rho_screw_pos_5
	rho_screw_pos_6 = rho_screw_pos_6
	rho_screw_pos_7 = rho_screw_pos_7
	rho_screw_pos_8 = rho_screw_pos_8
	rho_screw_pos_9 = rho_screw_pos_9
	rho_screw_pos_10 = rho_screw_pos_10
	rho_screw_pos_11 = rho_screw_pos_11
	rho_screw_pos_12 = rho_screw_pos_12
	rho_screw_neg_1 = rho_screw_neg_1
	rho_screw_neg_2 = rho_screw_neg_2
	rho_screw_neg_3 = rho_screw_neg_3
	rho_screw_neg_4 = rho_screw_neg_4
	rho_screw_neg_5 = rho_screw_neg_5
	rho_screw_neg_6 = rho_screw_neg_6
	rho_screw_neg_7 = rho_screw_neg_7
	rho_screw_neg_8 = rho_screw_neg_8
	rho_screw_neg_9 = rho_screw_neg_9
	rho_screw_neg_10 = rho_screw_neg_10
	rho_screw_neg_11 = rho_screw_neg_11
	rho_screw_neg_12 = rho_screw_neg_12
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
  end_time = 0.001 #0.1
  dt = 0.0005
  dtmin = 0.0001
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 1
  [../]
[]
