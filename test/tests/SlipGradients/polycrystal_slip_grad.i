# polycrystal
# additively manufactured 316L stainless steel parameters
# this model contains the slip gradients terms that
# induce growth of the GND density
# it is based on the paper by:
# Eralp Demir and Ivan Gutierrez-Urrutia
# Investigation of strain hardening near grain boundaries of an aluminum oligocrystal: 
# Experiments and crystal based finite element method
# International Journal of Plasticity
# Volume 136, January 2021, 102898

[GlobalParams]
  displacements = 'ux uy uz'
[]

[Mesh]
  [./polycrystal_generated_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 7
    nz = 1
    xmax = 19.0
    ymax = 13.3
    zmax = 1.9
    elem_type = HEX8
  [../]
  [./add_load_surface]
    type = SubdomainBoundingBoxGenerator
    input = polycrystal_generated_mesh
    bottom_left = '17.0 -0.001 -0.001'
	top_right = '19.1 14.0 2.0'
    block_id = 1
	block_name = 'elements_on_load_surface'
	location = 'INSIDE'
  [../]
[]

[Variables]
  [./ux]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  [./uy]
    order = FIRST
    family = LAGRANGE
  [../]

  [./uz]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz'
[]

[AuxVariables]

  # slip_increment is the slip rate
  # so units are 1/time
  [./slip_increment_1]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_2]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_3]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_4]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_5]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_6]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_7]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_8]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_9]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_10]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_11]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_12]
   order = FIRST
   family = MONOMIAL
  [../]

  [./slip_increment_vector]
    order = FIRST
    family = MONOMIAL
    components = 12
  [../]

  [./dslip_increment_dedge]
    order = CONSTANT
    family = MONOMIAL
    components = 12
  [../]
  
  [./dslip_increment_dscrew]
    order = CONSTANT
    family = MONOMIAL
    components = 12
  [../]
  
  [./rho_ssd_1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_ssd_12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_edge_12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./rho_gnd_screw_12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_12]
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

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'polycrystal.txt'
    nprop = 3
    read_type = element
  [../]
[]

[Functions]

  [./disp_load]
    type = ParsedFunction
    value = '19.0*0.01*t/60.0'
  [../]
  
  [./dts]
    type = PiecewiseConstant
    x = '0.0 1390.0'
    y = '0.1 0.1'
  [../]

[]

[AuxKernels]

  [./slip_increment_1]
   type = MaterialStdVectorAux
   variable = slip_increment_1
   property = slip_increment
   index = 0
   execute_on = timestep_end
  [../]
  [./slip_increment_2]
   type = MaterialStdVectorAux
   variable = slip_increment_2
   property = slip_increment
   index = 1
   execute_on = timestep_end
  [../]
  [./slip_increment_3]
   type = MaterialStdVectorAux
   variable = slip_increment_3   
   property = slip_increment
   index = 2
   execute_on = timestep_end
  [../]
  [./slip_increment_4]
   type = MaterialStdVectorAux
   variable = slip_increment_4
   property = slip_increment
   index = 3
   execute_on = timestep_end
  [../]
  [./slip_increment_5]
   type = MaterialStdVectorAux
   variable = slip_increment_5
   property = slip_increment
   index = 4
   execute_on = timestep_end
  [../]
  [./slip_increment_6]
   type = MaterialStdVectorAux
   variable = slip_increment_6
   property = slip_increment
   index = 5
   execute_on = timestep_end
  [../]
  [./slip_increment_7]
   type = MaterialStdVectorAux
   variable = slip_increment_7   
   property = slip_increment
   index = 6
   execute_on = timestep_end
  [../]
  [./slip_increment_8]
   type = MaterialStdVectorAux
   variable = slip_increment_8
   property = slip_increment
   index = 7
   execute_on = timestep_end
  [../]
  [./slip_increment_9]
   type = MaterialStdVectorAux
   variable = slip_increment_9
   property = slip_increment
   index = 8
   execute_on = timestep_end
  [../]
  [./slip_increment_10]
   type = MaterialStdVectorAux
   variable = slip_increment_10
   property = slip_increment
   index = 9
   execute_on = timestep_end
  [../]
  [./slip_increment_11]
   type = MaterialStdVectorAux
   variable = slip_increment_11   
   property = slip_increment
   index = 10
   execute_on = timestep_end
  [../]
  [./slip_increment_12]
   type = MaterialStdVectorAux
   variable = slip_increment_12
   property = slip_increment
   index = 11
   execute_on = timestep_end
  [../]

  [./build_slip_increment_vector]
    type = BuildArrayVariableAux
    variable = slip_increment_vector
    component_variables = 'slip_increment_1 slip_increment_2 slip_increment_3 slip_increment_4 slip_increment_5 slip_increment_6 slip_increment_7 slip_increment_8 slip_increment_9 slip_increment_10 slip_increment_11 slip_increment_12'
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
  
  [./rho_ssd_1]
    type = MaterialStdVectorAux
    variable = rho_ssd_1
    property = rho_ssd
    index = 0
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_2]
    type = MaterialStdVectorAux
    variable = rho_ssd_2
    property = rho_ssd
    index = 1
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_3]
    type = MaterialStdVectorAux
    variable = rho_ssd_3
    property = rho_ssd
    index = 2
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_4]
    type = MaterialStdVectorAux
    variable = rho_ssd_4
    property = rho_ssd
    index = 3
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_5]
    type = MaterialStdVectorAux
    variable = rho_ssd_5
    property = rho_ssd
    index = 4
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_6]
    type = MaterialStdVectorAux
    variable = rho_ssd_6
    property = rho_ssd
    index = 5
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_7]
    type = MaterialStdVectorAux
    variable = rho_ssd_7
    property = rho_ssd
    index = 6
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_8]
    type = MaterialStdVectorAux
    variable = rho_ssd_8
    property = rho_ssd
    index = 7
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_9]
    type = MaterialStdVectorAux
    variable = rho_ssd_9
    property = rho_ssd
    index = 8
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_10]
    type = MaterialStdVectorAux
    variable = rho_ssd_10
    property = rho_ssd
    index = 9
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_11]
    type = MaterialStdVectorAux
    variable = rho_ssd_11
    property = rho_ssd
    index = 10
    execute_on = timestep_end
  [../]
  
  [./rho_ssd_12]
    type = MaterialStdVectorAux
    variable = rho_ssd_12
    property = rho_ssd
    index = 11
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_1]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_1
    property = rho_gnd_edge
    index = 0
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_2]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_2
    property = rho_gnd_edge
    index = 1
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_3]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_3
    property = rho_gnd_edge
    index = 2
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_4]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_4
    property = rho_gnd_edge
    index = 3
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_5]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_5
    property = rho_gnd_edge
    index = 4
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_6]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_6
    property = rho_gnd_edge
    index = 5
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_7]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_7
    property = rho_gnd_edge
    index = 6
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_8]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_8
    property = rho_gnd_edge
    index = 7
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_9]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_9
    property = rho_gnd_edge
    index = 8
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_10]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_10
    property = rho_gnd_edge
    index = 9
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_11]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_11
    property = rho_gnd_edge
    index = 10
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_edge_12]
    type = MaterialStdVectorAux
    variable = rho_gnd_edge_12
    property = rho_gnd_edge
    index = 11
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_1]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_1
    property = rho_gnd_screw
    index = 0
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_2]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_2
    property = rho_gnd_screw
    index = 1
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_3]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_3
    property = rho_gnd_screw
    index = 2
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_4]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_4
    property = rho_gnd_screw
    index = 3
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_5]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_5
    property = rho_gnd_screw
    index = 4
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_6]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_6
    property = rho_gnd_screw
    index = 5
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_7]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_7
    property = rho_gnd_screw
    index = 6
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_8]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_8
    property = rho_gnd_screw
    index = 7
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_9]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_9
    property = rho_gnd_screw
    index = 8
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_10]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_10
    property = rho_gnd_screw
    index = 9
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_11]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_11
    property = rho_gnd_screw
    index = 10
    execute_on = timestep_end
  [../]
  
  [./rho_gnd_screw_12]
    type = MaterialStdVectorAux
    variable = rho_gnd_screw_12
    property = rho_gnd_screw
    index = 11
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_1]
    type = MaterialStdVectorAux
    variable = slip_resistance_1
    property = slip_resistance
    index = 0
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_2]
    type = MaterialStdVectorAux
    variable = slip_resistance_2
    property = slip_resistance
    index = 1
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_3]
    type = MaterialStdVectorAux
    variable = slip_resistance_3
    property = slip_resistance
    index = 2
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_4]
    type = MaterialStdVectorAux
    variable = slip_resistance_4
    property = slip_resistance
    index = 3
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_5]
    type = MaterialStdVectorAux
    variable = slip_resistance_5
    property = slip_resistance
    index = 4
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_6]
    type = MaterialStdVectorAux
    variable = slip_resistance_6
    property = slip_resistance
    index = 5
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_7]
    type = MaterialStdVectorAux
    variable = slip_resistance_7
    property = slip_resistance
    index = 6
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_8]
    type = MaterialStdVectorAux
    variable = slip_resistance_8
    property = slip_resistance
    index = 7
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_9]
    type = MaterialStdVectorAux
    variable = slip_resistance_9
    property = slip_resistance
    index = 8
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_10]
    type = MaterialStdVectorAux
    variable = slip_resistance_10
    property = slip_resistance
    index = 9
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_11]
    type = MaterialStdVectorAux
    variable = slip_resistance_11
    property = slip_resistance
    index = 10
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_12]
    type = MaterialStdVectorAux
    variable = slip_resistance_12
    property = slip_resistance
    index = 11
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
  [./z0_back]
    type = DirichletBC
    variable = uz
    boundary = back
    value = 0.0
  [../]

  [./y0_bottom]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0.0
  [../]
  
  [./x0_left]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0.0
  [../]

  [./x1_right]
    type = FunctionDirichletBC
    variable = ux
    boundary = right
    function = disp_load
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPGrain
    C_ijkl = '2.046e5 1.377e5 1.377e5 2.046e5 1.377e5 2.046e5 1.262e5 1.262e5 1.262e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
  [./stress]
    type = ComputeDislocationCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    maximum_substep_iteration = 2
	maxiter = 500
	maxiter_state_variable = 500
  [../]
  [./trial_xtalpl]
    type = CrystalPlasticityDislocationUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
	ao = 0.001
	xm = 0.1
	burgers_vector_mag = 0.000256
	shear_modulus = 86000.0 # MPa
	alpha_0 = 0.3
	r = 1.4
	tau_c_0 = 0.112
	k_0 = 0.02299282177563252
	y_c = 0.0019545318633428007 #0.0026 only edge annihilation
	init_rho_ssd = 35.925613042119906
	init_rho_gnd_edge = 0.0
	init_rho_gnd_screw = 0.0
	# These activate slip gradients
	# they are compulsory
	# codes currently has problems if not introduced
	# to remove the effect of slip gradients, zero arrays can be passed
	dslip_increment_dedge = dslip_increment_dedge
	dslip_increment_dscrew = dslip_increment_dscrew
	slip_increment_tolerance = 2.0
	stol = 0.1
	resistance_tol = 1.0
	print_state_variable_convergence_error_messages = true
  [../]
[]

[Postprocessors]
  [stress_on_load_surface]
    type = ElementAverageValue
    variable = stress_xx
	block = 'elements_on_load_surface'
    use_displaced_mesh = true
  []
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

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor -pc_hypre_boomeramg_print_statistics'
  petsc_options_value = 'hypre boomeramg 51 0.7 4 5 25 PMIS ext+i 2 0.3 0'

  line_search = 'none'

  automatic_scaling = true
  
  nl_abs_tol = 1e-5
  nl_rel_tol = 1e-5
  
  [./TimeStepper]
    type = FunctionDTGrowth
    function = dts
	cutback_factor_at_failure = 0.5
	growth_factor = 1.1
  [../]
  
  end_time = 0.1 #1338.0
  dtmin = 0.000001

[]

[Outputs]
  exodus = true
  csv = true
  interval = 1
[]
