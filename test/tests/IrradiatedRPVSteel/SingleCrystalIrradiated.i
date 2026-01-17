# model for irradiated RPV steel in:
# Ghiath Monnet, Ludovic Vincent, Lionel Gelebart
# Multiscale modeling of crystal plasticity in Reactor Pressure Vessel steels: Prediction of irradiation hardening
# Journal of Nuclear Materials 514 (2019) 128-138

[GlobalParams]
  displacements = 'ux uy uz'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
[]

[AuxVariables]
  [./fp_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho_ssd]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_resistance]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_increment]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./zero_first]
    order = FIRST
    family = MONOMIAL
  [../]
  [./slip_increment_vector]
    order = FIRST
    family = MONOMIAL
    components = 24
  [../]
  [./dslip_increment_dedge]
    order = CONSTANT
    family = MONOMIAL
    components = 24
  [../]
  [./dslip_increment_dscrew]
    order = CONSTANT
    family = MONOMIAL
    components = 24
  [../]
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz strain_xx strain_yy strain_zz strain_xy strain_xz strain_yz'
[]

[AuxKernels]
  [./fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = plastic_deformation_gradient
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = total_lagrangian_strain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  [./rho_ssd]
    type = MaterialStdVectorAux
    variable = rho_ssd
    property = rho_ssd
    index = 4
    execute_on = timestep_end
  [../]
  [./slip_resistance]
    type = MaterialStdVectorAux
	variable = slip_resistance
	property = slip_resistance
	index = 4
	execute_on = timestep_end
  [../]
  [./slip_increment]
    type = MaterialStdVectorAux
	variable = slip_increment
	property = slip_increment
	index = 4
	execute_on = timestep_end
  [../]
  
  # slip gradients are set to zero
  [./build_slip_increment_vector]
    type = BuildArrayVariableAux
    variable = slip_increment_vector
    component_variables = 'zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first zero_first'
  [../]
  [./edge_directional_derivative]
    type = ArrayDirectionalDerivative
    variable = dslip_increment_dedge
    gradient_variable = slip_increment_vector
    dislo_character = edge
  	execute_on = timestep_end
  [../]
  [./screw_directional_derivative]
    type = ArrayDirectionalDerivative
    variable = dslip_increment_dscrew
    gradient_variable = slip_increment_vector
    dislo_character = screw
  	execute_on = timestep_end
  [../]  
[]

[BCs]
  [./symmy]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0
  [../]
  [./symmx]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0
  [../]
  [./symmz]
    type = DirichletBC
    variable = uz
    boundary = back
    value = 0
  [../]
  [./tdisp]
    type = FunctionDirichletBC
    variable = uz
    boundary = front
    function = '0.001*t'
  [../]
[]

[UserObjects]
  
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'EulerAngles.txt'
    nprop = 3
    read_type = element
  [../]
  
[]

[Materials]

  # Elastic constants from Table 2 in:
  # Sudook A.Kim, Ward L.Johnson
  # Elastic constants and internal friction of martensitic steel, ferritic-pearlitic steel, and α-iron
  # Materials Science and Engineering: A
  # Volumes 452–453, 15 April 2007, Pages 633-639
  # https://www.sciencedirect.com/science/article/pii/S0921509306026244
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPGrain
    C_ijkl = '2.3e5 1.35e5 1.35e5 2.3e5 1.35e5 2.3e5 1.17e5 1.17e5 1.17e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
  
  [./stress]
    type = ComputeDislocationCrystalPlasticityStress
    crystal_plasticity_models = 'irradiatedRPVsteel'
    tan_mod_type = exact
    maximum_substep_iteration = 5
	rtol = 1e-4
	abs_tol = 1e-4
  [../]
  
  # Parameters are like in the original article
  # but d_grain = 50.0 is recalibrated to match the 464 MPa
  # L. Vincent, M. Libert, B. Marini, C. Reyb  
  # Towards a modelling of RPV steel brittle fracture using crystal plasticity
  # computations on polycrystalline aggregates
  # Journal of Nuclear Materials, Volume 406, Issue 1, 1 November 2010, Pages 91-96
  
  [./irradiatedRPVsteel]
    type = CrystalPlasticityIrradiatedRPVSteel
    number_slip_systems = 24
    slip_sys_file_name = 'input_slip_sys_bcc24.txt'
    resistance_tol = 0.1
	slip_increment_tolerance = 0.1
	burgers_vector_mag = 0.000256
	shear_modulus = 86000.0
	RT_shear_modulus = 86000.0
	a_self = 0.1
	a_col = 0.7
	K_Hall_Petch = 480.0
	d_grain = 50.0 #6.9 micron
	rho_carbide = 0.0608
	a_carbide = 0.0
	C_DL_diameter = 0.0256
	a_DL = 0.25
	C_SC_diameter = 0.0256
	a_SC = 0.04
	rho_ref = 1.0
	ao = 0.000001
	xm = 0.01
	attack_frequency = 2.0e11
	minimum_screw_length = 0.010
	Gibbs_free_energy_slip = 0.84
	k = 8.6e-5
	const_slip_resistance_110 = 360.0
	const_slip_resistance_112_TW = 410.0
	const_slip_resistance_112_AT = 480.0
	K_self = 17.0
	K_forest = 5.666
	y_drag = 0.002
	lambda_DL = 1.0
	lambda_SC = 1.0
	init_rho_ssd = 10.0
	init_rho_gnd_edge = 0.0
	init_rho_gnd_screw = 0.0
	
	# To match the yield strength of 540 MPa of the irradiated Euromaterial A, 
	# assuming the same initial density of solute clusters and dislocation loops,
	# the following are the resulting initial densities of dislocation loops due to irradiation
	# and solute clusters due to irradiation
	
	init_C_DL = 2600.0 # 1/micron^3
	init_C_SC = 2600.0 # 1/micron^3
	
	rho_tol = 1.0
	dslip_increment_dedge = dslip_increment_dedge
	dslip_increment_dscrew = dslip_increment_dscrew
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
  nl_abs_tol = 1e-6
  nl_rel_step_tol = 1e-6
  dtmax = 1.0
  nl_rel_tol = 1e-6
  dtmin = 0.000001
  dt = 0.01
  end_time = 0.01 # end_time = 50.0 to reach 5% strain
  nl_abs_step_tol = 1e-6
[]

[Outputs]
  exodus = true
  time_step_interval = 1 # 10
[]
