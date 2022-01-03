# This is the simulation used in
# Nicolò Grilli, Daijun Hu, Dewen Yushu, Fan Chen, Wentao Yan
# Crystal plasticity model of residual stress in additive manufacturing 
# using the element elimination and reactivation method
# Computational Mechanics (Springer) 2021
# https://doi.org/10.1007/s00466-021-02116-z
# https://link.springer.com/article/10.1007/s00466-021-02116-z
# it is used in Figures 8, 11, 12, 14 of the article

# the temperature file Temperature4umToRT2.txt
# is truncated to two temperature time steps
# so that it does not occupy much space in github

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 50 
    ny = 40
    nz = 27
    xmax = 200.0
    ymax = 160.0
    zmax = 108.0
    elem_type = HEX8
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./active_domain]
    input = gen
    type = SubdomainBoundingBoxGenerator
    bottom_left = '0.0 0.0 0.0'
    top_right = '200.0 160.0 104.0'
    block_id = 1
  [../]
  [./deactivated_domain]
    input = active_domain
    type = SubdomainBoundingBoxGenerator
    bottom_left = '0.0 0.0 104.0'
    top_right = '200.0 160.0 108.0'
    block_id = 2
  [../]
  [./sidesets]
    input = deactivated_domain
    type = SideSetsAroundSubdomainGenerator
    normal = '0 0 1'
    block = 1
    new_boundary = 'moving_interface'
  []
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
    [./InitialCondition]
      type = FunctionIC
      function = disp_z_init
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
  
  [./phase_temp]
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
    block = '1'
  [../]

  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./stress_yx]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./stress_zx]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./stress_zy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./e_xx]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./fp_xx]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./fp_xy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./fp_xz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./fp_yx]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./fp_yz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./fp_zx]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./fp_zy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./fp_zz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./euler1]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./euler2]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./euler3]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./lattice_strain_xx]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./lattice_strain_xy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./lattice_strain_xz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./lattice_strain_yx]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./lattice_strain_yy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./lattice_strain_yz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./lattice_strain_zx]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./lattice_strain_zy]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./lattice_strain_zz]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss1]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss2]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss3]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss4]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss5]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
 
  [./gss6]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss7]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss8]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss9]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss10]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss11]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  
  [./gss12]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
[]

[Functions]
  [./dts]
    type = PiecewiseConstant
    x = '0   44.0 44.5 45.5 46.0 89.0 89.5 90.5 91.0 134.0 134.5 135.5 136.0 179.0 179.5 180.5 181.0 224.0 224.5 225.5 226.0 269.0 269.5 270.5 271.0 314.0 314.5 315.5 316.0 359.0 359.5 360.5 361.0 404.0 404.5 405.5 406.0 449.0 449.5 450.5 451.0 494.0 494.5 495.5 496.0 539.0 539.5 540.5 541.0 584.0 584.5 585.5 586.0 629.0 629.5 630.5 631.0 674.0 674.5 675.5 676.0 719.0 719.5 720.5 721.0 764.0 764.5 765.5 766.0 809.0 809.5 810.5 811.0 854.0 854.5 855.5 856.0 899.0 899.5 900.5 901.0 3455.0 3505.0 3555.0'
    y = '1.0  0.1  0.1  0.1  1.0  0.1  0.1  0.1  1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1  10.0   10.0    1.0    1.0'
  [../]
  [./disp_z_init]
    type = ParsedFunction
    value = 'if(t/20.0,3.0,0.0)'	
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

  [./temperature_read]
    type = LaserTempReadFile
    temperature_file_name = 'Temperature4umToRT2.txt'
    temperature_num_step = 2
  [../]
  
  [./activated_elem_uo]
    type = ActDeactElementsMelting
    execute_on = timestep_begin
    temperature = phase_temp
    active_subdomain_id = 1
	deactive_subdomain_id = 2
    expand_boundary_name = 'moving_interface'
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = false
    add_variables = true
    block = '1'
  [../]
  [./time_derivative_x]
    type = TimeDerivative
    variable = disp_x
	block = '2'
  [../]
  [./time_derivative_y]
    type = TimeDerivative
    variable = disp_y
	block = '2'
  [../]
  [./time_derivative_z]
    type = TimeDerivative
    variable = disp_z
	block = '2'
  [../]
  [./reaction_z]
    type = Reaction
    variable = disp_z
	rate = 0.2
	block = '2'
  [../]
  [./force_z]
    type = BodyForce
    variable = disp_z
	value = 0.6
	block = '2'
  [../]
[]

[AuxKernels]

  [./tempfuncaux]
    type = LaserTempReadFileAux
    variable = temp
    temperature_read_user_object = temperature_read
    temperature_time_step = 45.0
	degrade_eigenstrain = true
  [../]
  
  [./phasetempfuncaux]
    type = TempActDeactElemsCFDStepsAux
    variable = phase_temp
    temperature_read_user_object = temperature_read
    temperature_time_step = 45.0
  [../]

  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
    index_j = 1
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./stress_xz]
    type = RankTwoAux
    variable = stress_xz
    rank_two_tensor = stress
    index_j = 2
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./stress_yx]
    type = RankTwoAux
    variable = stress_yx
    rank_two_tensor = stress
    index_j = 0
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./stress_yz]
    type = RankTwoAux
    variable = stress_yz
    rank_two_tensor = stress
    index_j = 2
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]
 
  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./stress_zx]
    type = RankTwoAux
    variable = stress_zx
    rank_two_tensor = stress
    index_j = 0
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./stress_zy]
    type = RankTwoAux
    variable = stress_zy
    rank_two_tensor = stress
    index_j = 1
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]

  [./gss1]
    type = MaterialStdVectorAux
    variable = gss1
    property = gss
    index = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./gss2]
    type = MaterialStdVectorAux
    variable = gss2
    property = gss
    index = 1
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./gss3]
    type = MaterialStdVectorAux
    variable = gss3
    property = gss
    index = 2
    execute_on = timestep_end
    block = '1'
  [../]
 
  [./gss4]
    type = MaterialStdVectorAux
    variable = gss4
    property = gss
    index = 3
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./gss5]
    type = MaterialStdVectorAux
    variable = gss5
    property = gss
    index = 4
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./gss6]
    type = MaterialStdVectorAux
    variable = gss6
    property = gss
    index = 5
    execute_on = timestep_end
    block = '1'
  [../]
 
  [./gss7]
    type = MaterialStdVectorAux
    variable = gss7
    property = gss
    index = 6
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./gss8]
    type = MaterialStdVectorAux
    variable = gss8
    property = gss
    index = 7
    execute_on = timestep_end
    block = '1'
  [../]
 
  [./gss9]
    type = MaterialStdVectorAux
    variable = gss9
    property = gss
    index = 8
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./gss10]
    type = MaterialStdVectorAux
    variable = gss10
    property = gss
    index = 9
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./gss11]
    type = MaterialStdVectorAux
    variable = gss11
    property = gss
    index = 10
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./gss12]
    type = MaterialStdVectorAux
    variable = gss12
    property = gss
    index = 11
    execute_on = timestep_end
    block = '1'
  [../]

  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = lage
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./e_xx]
    type = RankTwoAux
    variable = e_xx
    rank_two_tensor = lage
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = lage
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]

  [./fp_xx]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = fp
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./fp_xy]
    type = RankTwoAux
    variable = fp_xy
    rank_two_tensor = fp
    index_j = 1
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./fp_xz]
    type = RankTwoAux
    variable = fp_xz
    rank_two_tensor = fp
    index_j = 2
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./fp_yx]
    type = RankTwoAux
    variable = fp_yx
    rank_two_tensor = fp
    index_j = 0
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]

  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]

  [./fp_yz]
    type = RankTwoAux
    variable = fp_yz
    rank_two_tensor = fp
    index_j = 2
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]

  [./fp_zx]
    type = RankTwoAux
    variable = fp_zx
    rank_two_tensor = fp
    index_j = 0
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]

  [./fp_zy]
    type = RankTwoAux
    variable = fp_zy
    rank_two_tensor = fp
    index_j = 1
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]

  [./fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = fp
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]

  [./euler1]
    type = MaterialRealVectorValueAux
    variable = euler1
    property = Euler_angles
    component = 0
    execute_on = timestep_end
    block = '1'
  [../]

  [./euler2]
    type = MaterialRealVectorValueAux
    variable = euler2
    property = Euler_angles
    component = 1
    execute_on = timestep_end
    block = '1'
  [../]

  [./euler3]
    type = MaterialRealVectorValueAux
    variable = euler3
    property = Euler_angles
    component = 2
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./lattice_strain_xx]
    type = RankTwoAux
    variable = lattice_strain_xx
    rank_two_tensor = lattice_strain
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]

  [./lattice_strain_zz]
   type = RankTwoAux
    variable = lattice_strain_zz
    rank_two_tensor = lattice_strain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]

  [./lattice_strain_yy]
   type = RankTwoAux
    variable = lattice_strain_yy
    rank_two_tensor = lattice_strain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./lattice_strain_xy]
   type = RankTwoAux
    variable = lattice_strain_xy
    rank_two_tensor = lattice_strain
    index_j = 1
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./lattice_strain_yx]
   type = RankTwoAux
    variable = lattice_strain_yx
    rank_two_tensor = lattice_strain
    index_j = 0
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./lattice_strain_yz]
   type = RankTwoAux
    variable = lattice_strain_yz
    rank_two_tensor = lattice_strain
    index_j = 2
    index_i = 1
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./lattice_strain_zy]
   type = RankTwoAux
    variable = lattice_strain_zy
    rank_two_tensor = lattice_strain
    index_j = 1
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./lattice_strain_xz]
   type = RankTwoAux
    variable = lattice_strain_xz
    rank_two_tensor = lattice_strain
    index_j = 2
    index_i = 0
    execute_on = timestep_end
    block = '1'
  [../]
  
  [./lattice_strain_zx]
   type = RankTwoAux
    variable = lattice_strain_zx
    rank_two_tensor = lattice_strain
    index_j = 0
    index_i = 2
    execute_on = timestep_end
    block = '1'
  [../]
[]

[BCs]
  [./z0_back]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
	use_displaced_mesh = false
  [../]

  [./y0]
    type = DirichletBC
    variable = disp_y
    boundary = 'back bottom top left right'
    value = 0.0
	use_displaced_mesh = false
  [../]

  [./x0]
    type = DirichletBC
    variable = disp_x
    boundary = 'back bottom top left right'
    value = 0.0
	use_displaced_mesh = false
  [../]
[]

[Postprocessors]

[]

[Materials]
  [./crysp]
    type = FiniteStrainCrystalPlasticityThermal
    block = '1'
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 12 #Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 12 0.0000000000001 0.1' # slip rate equations parameters
    hprops = '1.0 3839.0 213.0 302.0 2.5' # hardening properties
    gprops = '1 12 213.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = none
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
# Calibrated using table 1 in:
# W.Jiang, Y.Zhang and W.Woo
# Using heat sink technology to decrease residual stress
# in 316L stainless steel welding joint:Finite element simulation
# Int.J.Press.Vessel.Pip.
# VOLUME 92, pp.56-62, 2012
    dCTE_dT = '0.0' #'0.01011e-6'
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
    residual_stiffness = 0.01
    temperature_read_user_object = temperature_read
    temperature_time_step = 45.0
    reference_temperature = '303.0'
    activate_elems = true
    block = '1'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
    block = '1'
  [../]
  [./dummy_mat_inactive]
    type = GenericConstantMaterial
    prop_names = 'dummy_mat'
    prop_values = '0.0'
  	block = '2'
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
  #petsc_options = '-snes_ksp_ew -snes_test_jacobian'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  line_search = 'default'
  automatic_scaling = true

  l_max_its = 30
  nl_max_its = 10
  nl_rel_tol = 1e-3
  nl_abs_tol = 1e-3
  
  l_tol = 1e-5

  [./TimeStepper]
    type = FunctionDTGrowth
    function = dts
	cutback_factor_at_failure = 0.1
	growth_factor = 1.5
  [../]
  
  start_time = 0.0
  end_time = 1.0 # 3555.0 microseconds

  dtmin = 1.0e-99
  timestep_tolerance = 1.0e-99  

[]

[Outputs]
  csv = false
  #[./my]
  #  type = Checkpoint
  #  num_files = 2
  #  interval = 45
  #[../]
  [./out]
    type = Exodus
    interval = 1
  [../]
[]

