[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  xmax = 2.0
  ymax = 2.0
  zmax = 2.0
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
    # Use the initial Condition block underneath the variable
    # for which we want to apply this initial condition
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
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
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = ConstantIC
      value = '298.0'
    [../]
  [../]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./euler1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./euler2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./euler3]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]

[]

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'GGEuler.txt'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    read_type = element
  [../]
  [./temperature_read]
    type = LaserTempReadFile
	temperature_file_name = 'temperature.txt'
	temperature_num_step = 3
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
  [./tempfuncaux]
    type = LaserTempReadFileAux
    variable = temp
    temperature_read_user_object = temperature_read
	temperature_time_step = 0.1
    block = 'ANY_BLOCK_ID 0'
  [../]

  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]

  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]  

  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]

  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  [../]

  [./stress_xz]
    type = RankTwoAux
    variable = stress_xz
    rank_two_tensor = stress
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  [../]

  [./stress_yz]
    type = RankTwoAux
    variable = stress_yz
    rank_two_tensor = stress
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  [../]

  [./euler1]
    type = MaterialRealVectorValueAux
    variable = euler1
    property = Euler_angles
    component = 0
    execute_on = timestep_end
  [../]
  
  [./euler2]
    type = MaterialRealVectorValueAux
    variable = euler2
    property = Euler_angles
    component = 1
    execute_on = timestep_end
  [../]
  
  [./euler3]
    type = MaterialRealVectorValueAux
    variable = euler3
    property = Euler_angles
    component = 2
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./z0_bot]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./y0_bot]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./x0_bot]
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
    gtol = 1e-2
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 12 #Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' # slip rate equations parameters 
    hprops = '1.0 3839.0 213.0 302.0 2.5' # hardening properties
    gprops = '1 12 213.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    thermal_expansion = '44.73e-6'
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
# Calibrated using table 1 in:
# W.Jiang, Y.Zhang and W.Woo
# Using heat sink technology to decrease residual stress
# in 316L stainless steel welding joint:Finite element simulation
# Int.J.Press.Vessel.Pip.
# VOLUME 92, pp.56-62, 2012
    dCTE_dT='0.01011e-6'
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
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-5
  l_tol = 1e-5
 
  start_time = 0.0
  end_time = 0.2
  dt = 0.05
  dtmin = 0.001

[]

[Outputs]
  csv = false
  [./out]
    type = Exodus
    interval = 1
  [../]
[]
