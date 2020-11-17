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
[]

[AuxVariables]
  [./temp]
    order = FIRST
    family = LAGRANGE
	[./InitialCondition]
      type = ConstantIC
      value = 1000.0
    [../]
  [../]

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./e_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./gss]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]

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
	temperature_num_step = 8
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    add_variables = true
  [../]
[]

[AuxKernels]

  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]  

  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

 [./fp_xx]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = fp
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

 [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = fp
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = lage
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./e_xx]
    type = RankTwoAux
    variable = e_xx
    rank_two_tensor = lage
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = lage
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./tempfuncaux]
    type = LaserTempReadFileAux
    variable = temp
    temperature_read_user_object = temperature_read
	temperature_time_step = 1.0
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./gss]
    type = MaterialStdVectorAux
    variable = gss
    property = gss
    index = 0
    execute_on = timestep_end
    block = 'ANY_BLOCK_ID 0'
  [../]
[]

[BCs]

  [./z_bot]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./y_bot]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./x_bot]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

[]
 
[Postprocessors]

[]

[Materials]
  [./crysp]
    type = FiniteStrainCrystalPlasticityThermal
    block = 0
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
    hprops = '1.0 541.5 216.0 250.0 2.5' # hardening properties
    gprops = '1 12 216.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    thermal_expansion = '17e-6'
    reference_temperature = '1000.0'
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
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorMelting
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
	residual_stiffness = 0.001
	temperature_read_user_object = temperature_read
	temperature_time_step = 1.0
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
  end_time = 0.02
  dt = 0.01
  dtmin = 0.01
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 1
  [../]
[]
