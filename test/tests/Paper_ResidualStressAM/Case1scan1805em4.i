# Case 1 simulation used in figure 10 in 
# Daijun Hu, Nicolò Grilli, Lu Wang, Min Yang, Wentao Yan
# Microscale residual stresses in additively manufactured stainless steel: 
# Computational simulation
# Journal of the Mechanics and Physics of Solids
# Volume 161, April 2022, 104822
# Temperature file has been reduced to 2 temperature time steps
# to reduce the size on github account
# This file is for the laser scan simulation

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 48
  ny = 53
  nz = 38
  xmax = 192.0
  ymax = 212.0
  zmax = 152.0
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
  [./vel_x]
  [../]
  [./accel_x]
  [../]
  [./vel_y]
  [../]
  [./accel_y]
  [../]
  [./vel_z]
  [../]
  [./accel_z]
  [../]
  
  [./temp]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = ConstantIC
      value = '303.0'
    [../]
  [../]

  [./stress_xx]
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

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zy]
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
  [./fp_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_yx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_zx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_zy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_zz]
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

  [./lattice_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./lattice_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./lattice_strain_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./lattice_strain_yx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./lattice_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./lattice_strain_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./lattice_strain_zx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./lattice_strain_zy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./lattice_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss12]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./dts]
    type = PiecewiseLinear
    x = '0.00000 0.03200 0.07920'
    y = '0.000035 0.000035 0.00020'
  [../]
[]

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'GGEuler4umChen2019.txt'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    read_type = element
  [../]

  # reduced to 2 temperature steps
  # but for the full simulation 198 temperature steps were used   
  [./temperature_read]
    type = LaserTempReadFile
    temperature_file_name = 'Temperature4umToRT2.txt'
    temperature_num_step = 2 #198
  [../]
[]

[Kernels]
  [./DynamicTensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = false
    add_variables = true
    alpha = -0.3
    zeta = 0.0001
  [../]
  [./inertia_x]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = 0.3025
    gamma = 0.6
    alpha = -0.3
    eta = 0.0001
  [../]
  [./inertia_y]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = 0.3025
    gamma = 0.6
    alpha = -0.3
    eta = 0.0001
  [../]
  [./inertia_z]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_z
    velocity = vel_z
    acceleration = accel_z
    beta = 0.3025
    gamma = 0.6
    alpha = -0.3
    eta = 0.0001
  [../]
[]

[AuxKernels]

  [./accel_x]
    type = NewmarkAccelAux
	use_displaced_mesh = false
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = 0.3025
    execute_on = timestep_end
  [../]
  [./vel_x]
    type = NewmarkVelAux
	use_displaced_mesh = false
    variable = vel_x
    acceleration = accel_x
    gamma = 0.6
    execute_on = timestep_end
  [../]
  [./accel_y]
    type = NewmarkAccelAux
	use_displaced_mesh = false
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = 0.3025
    execute_on = timestep_end
  [../]
  [./vel_y]
    type = NewmarkVelAux
	use_displaced_mesh = false
    variable = vel_y
    acceleration = accel_y
    gamma = 0.6
    execute_on = timestep_end
  [../]
  [./accel_z]
    type = NewmarkAccelAux
	use_displaced_mesh = false
    variable = accel_z
    displacement = disp_z
    velocity = vel_z
    beta = 0.3025
    execute_on = timestep_end
  [../]
  [./vel_z]
    type = NewmarkVelAux
	use_displaced_mesh = false
    variable = vel_z
    acceleration = accel_z
    gamma = 0.6
    execute_on = timestep_end
  [../]
  
  # each temperature time step corresponds to 400 microseconds
  [./tempfuncaux]
    type = LaserTempReadFileAux
    variable = temp
    temperature_read_user_object = temperature_read
    temperature_time_step = 0.000400
  [../]

  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
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
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./stress_yx]
    type = RankTwoAux
    variable = stress_yx
    rank_two_tensor = stress
    index_j = 0
    index_i = 1
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
  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  [./stress_zx]
    type = RankTwoAux
    variable = stress_zx
    rank_two_tensor = stress
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  [../]
  [./stress_zy]
    type = RankTwoAux
    variable = stress_zy
    rank_two_tensor = stress
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  [../]

  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = lage
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  [./e_xx]
    type = RankTwoAux
    variable = e_xx
    rank_two_tensor = lage
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = lage
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]

  [./fp_xx]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = fp
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./fp_xy]
    type = RankTwoAux
    variable = fp_xy
    rank_two_tensor = fp
    index_j = 1
    index_i = 0
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
  [./fp_yx]
    type = RankTwoAux
    variable = fp_yx
    rank_two_tensor = fp
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  [../]

  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]

  [./fp_yz]
    type = RankTwoAux
    variable = fp_yz
    rank_two_tensor = fp
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  [../]

  [./fp_zx]
    type = RankTwoAux
    variable = fp_zx
    rank_two_tensor = fp
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  [../]

  [./fp_zy]
    type = RankTwoAux
    variable = fp_zy
    rank_two_tensor = fp
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  [../]

  [./fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = fp
    index_j = 2
    index_i = 2
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
  [./lattice_strain_xx]
    type = RankTwoAux
    variable = lattice_strain_xx
    rank_two_tensor = lattice_strain
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]

  [./lattice_strain_zz]
   type = RankTwoAux
    variable = lattice_strain_zz
    rank_two_tensor = lattice_strain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]

  [./lattice_strain_yy]
   type = RankTwoAux
    variable = lattice_strain_yy
    rank_two_tensor = lattice_strain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./lattice_strain_xy]
   type = RankTwoAux
    variable = lattice_strain_xy
    rank_two_tensor = lattice_strain
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  [../]
  [./lattice_strain_yx]
   type = RankTwoAux
    variable = lattice_strain_yx
    rank_two_tensor = lattice_strain
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  [../]
  [./lattice_strain_yz]
   type = RankTwoAux
    variable = lattice_strain_yz
    rank_two_tensor = lattice_strain
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  [../]
  [./lattice_strain_zy]
   type = RankTwoAux
    variable = lattice_strain_zy
    rank_two_tensor = lattice_strain
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  [../]
  [./lattice_strain_xz]
   type = RankTwoAux
    variable = lattice_strain_xz
    rank_two_tensor = lattice_strain
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  [../]
  [./lattice_strain_zx]
   type = RankTwoAux
    variable = lattice_strain_zx
    rank_two_tensor = lattice_strain
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  [../]

  [./gss1]
    type = MaterialStdVectorAux
    variable = gss1
    property = gss
    index = 0
    execute_on = timestep_end
  [../]
  [./gss2]
    type = MaterialStdVectorAux
    variable = gss2
    property = gss
    index = 1
    execute_on = timestep_end
  [../]
  [./gss3]
    type = MaterialStdVectorAux
    variable = gss3
    property = gss
    index = 2
    execute_on = timestep_end
  [../]
  [./gss4]
    type = MaterialStdVectorAux
    variable = gss4
    property = gss
    index = 3
    execute_on = timestep_end
  [../]
  [./gss5]
    type = MaterialStdVectorAux
    variable = gss5
    property = gss
    index = 4
    execute_on = timestep_end
  [../]
  [./gss6]
    type = MaterialStdVectorAux
    variable = gss6
    property = gss
    index = 5
    execute_on = timestep_end
  [../]
  [./gss7]
    type = MaterialStdVectorAux
    variable = gss7
    property = gss
    index = 6
    execute_on = timestep_end
  [../]
  [./gss8]
    type = MaterialStdVectorAux
    variable = gss8
    property = gss
    index = 7
    execute_on = timestep_end
  [../]
  [./gss9]
    type = MaterialStdVectorAux
    variable = gss9
    property = gss
    index = 8
    execute_on = timestep_end
  [../]
  [./gss10]
    type = MaterialStdVectorAux
    variable = gss10
    property = gss
    index = 9
    execute_on = timestep_end
  [../]
  [./gss11]
    type = MaterialStdVectorAux
    variable = gss11
    property = gss
    index = 10
    execute_on = timestep_end
  [../]
  [./gss12]
    type = MaterialStdVectorAux
    variable = gss12
    property = gss
    index = 11
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

  [./y0_bott]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./y0_top]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 0.0
  [../]

  [./x0_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  
  [./x0_right]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0.0
  [../]

[]

[Postprocessors]

[]

[Materials]
  [./density]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'density'
    prop_values = '7.87e-15'
  [../]
  [./crysp]
    type = FiniteStrainCrystalPlasticityThermal
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 12 #Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 12 0.0005 0.1' # slip rate equations parameters
    hprops = '1.0 3839.0 180.0 302.0 2.5' # hardening properties
    gprops = '1 12 180.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    thermal_expansion = '44.73e-6'
    reference_temperature = '303.0'
    temp = temp
    slip_incr_tol = 0.001
    maxiter = 500
    maxitergss = 500
    maximum_substep_iteration = 3
    gen_random_stress_flag = true
    rtol = '1.0e-4'
    abs_tol = '1.0e-4'
# Calibrated using table 1 in:
# M.R. DAYMOND and P.J. BOUCHARD
# Elastoplastic Deformation of 316 Stainless Steel Under
# Tensile Loading at Elevated Temperatures
# METALLURGICAL AND MATERIALS TRANSACTIONS A
# VOLUME 37A, JUNE 2006â€?873
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
    type = ComputeElasticityTensorMelting
# Elastic constants of 316L SS from:
# Clausen, B., Lorentzen, T. and Leffers, T.
# Self-consistent modelling of the plastic
# deformation of FCC polycrystals and its implications for diffraction
# measurements of internal stresses.
# Acta Mater. 46, 3087â€?098 (1998).
    C_ijkl = '2.046e5 1.377e5 1.377e5 2.046e5 1.377e5 2.046e5 1.262e5 1.262e5 1.262e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
    temp = temp
    dC11_dT = 0.0004415
    dC12_dT = 0.0003275
    dC44_dT = 0.0004103
        residual_stiffness = 0.01
        temperature_read_user_object = temperature_read
        temperature_time_step = 0.000400
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
  petsc_options_value = 'hypre    boomeramg          51'
  line_search = 'none'
  
  l_max_its = 17
  nl_max_its = 15
  
  nl_rel_tol = 1e-3
  nl_abs_tol = 1e-3
  
  #l_tol = 2e-6

  [./TimeStepper]
    type = FunctionDTGrowth
    function = dts
	cutback_factor_at_failure = 0.1
	growth_factor = 1.2  
  [../]
  
  # 78 ms is the total time of the full simulation
  start_time = 0.0
  end_time = 0.000035 #0.078
  dtmin = 1.0e-20

[]

[Outputs]
  csv = false
  # checkpoints were saved every 20 steps
  # in the complete simulation
  # restart is needed for the stress relaxation simulation
  [./my]
    type = Checkpoint
    num_files = 5
    interval = 1 #20
  [../]
  [./out]
    type = Exodus
    interval = 1 #55
    checkpoint = true
  [../]
[]

