# simplified bcc model for irradiated RPV steel
# master app
# neper generated geometry with hexahedra
# JRQ steel model parameters
# crystal plasticity solved with undamaged stress to improve convergence

[GlobalParams]
  displacements = 'ux uy uz'
[]

# neper generated mesh with 8 grains
[Mesh]
  [./nepermesh]
    type = FileMeshGenerator
    file = 'n8-id1.msh'
  [../]
  [./left_modifier]
    type = BoundingBoxNodeSetGenerator
    input = nepermesh
    new_boundary = left
    top_right = '0.001 40.001 1.001'
    bottom_left = '-0.001 -0.001 -0.001'
  [../]
  [./bottom_modifier]
    type = BoundingBoxNodeSetGenerator
    input = left_modifier
    new_boundary = bottom
    top_right = '40.001 0.001 1.001'
    bottom_left = '-0.001 -0.001 -0.001'
  [../]
  [./back_modifier]
    type = BoundingBoxNodeSetGenerator
    input = bottom_modifier
    new_boundary = back
    top_right = '40.001 40.001 0.001'
    bottom_left = '-0.001 -0.001 -0.001'
  [../]
  [./right_modifier]
    type = BoundingBoxNodeSetGenerator
    input = back_modifier
    new_boundary = right
    top_right = '40.001 40.001 1.001'
    bottom_left = '39.999 -0.001 -0.001'
  [../]
  [./front_modifier]
    type = BoundingBoxNodeSetGenerator
    input = right_modifier
    new_boundary = front
    top_right = '40.001 40.001 1.001'
    bottom_left = '-0.001 -0.001 0.999'
  [../]
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    execute_on = 'TIMESTEP_END'
  []
[]

[AuxVariables]

  [./euler1]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./slip_rate_1]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_2]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_3]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_4]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_5]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_6]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_7]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_8]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_9]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_10]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_11]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_12]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_13]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_14]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_15]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_16]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_17]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_18]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_19]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_20]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_21]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_22]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_23]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_rate_24]
   order = FIRST
   family = MONOMIAL
  [../]  
  
  # build slip rate vector
  [./slip_rate_vector]
    order = FIRST
    family = MONOMIAL
    components = 24
  [../]
  
  # derivatives of the slip rate
  # with respect to the edge and screw
  # direction of motion
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
  
  [./c]
  [../]
  
  [./elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./delastic_energy_dc]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./d2elastic_energy_dc2]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz'
[]

# no off diagonal Jacobian kernels
# because mechanics and damage problems are independent

[AuxKernels]
 
  [./euler1]
    type = MaterialRealVectorValueAux
    variable = euler1
    property = Euler_angles
    component = 0
    execute_on = timestep_end
  [../]

  # slip rate
  # of each slip system
  # so units are 1/time
  # the corresponding material property in the material class
  # is called slip_increment but it corresponds to the rate
  [./slip_rate_1]
   type = MaterialStdVectorAux
   variable = slip_rate_1
   property = slip_increment
   index = 0
   execute_on = timestep_end
  [../]
  [./slip_rate_2]
   type = MaterialStdVectorAux
   variable = slip_rate_2
   property = slip_increment
   index = 1
   execute_on = timestep_end
  [../]
  [./slip_rate_3]
   type = MaterialStdVectorAux
   variable = slip_rate_3   
   property = slip_increment
   index = 2
   execute_on = timestep_end
  [../]
  [./slip_rate_4]
   type = MaterialStdVectorAux
   variable = slip_rate_4
   property = slip_increment
   index = 3
   execute_on = timestep_end
  [../]
  [./slip_rate_5]
   type = MaterialStdVectorAux
   variable = slip_rate_5
   property = slip_increment
   index = 4
   execute_on = timestep_end
  [../]
  [./slip_rate_6]
   type = MaterialStdVectorAux
   variable = slip_rate_6
   property = slip_increment
   index = 5
   execute_on = timestep_end
  [../]
  [./slip_rate_7]
   type = MaterialStdVectorAux
   variable = slip_rate_7   
   property = slip_increment
   index = 6
   execute_on = timestep_end
  [../]
  [./slip_rate_8]
   type = MaterialStdVectorAux
   variable = slip_rate_8
   property = slip_increment
   index = 7
   execute_on = timestep_end
  [../]
  [./slip_rate_9]
   type = MaterialStdVectorAux
   variable = slip_rate_9
   property = slip_increment
   index = 8
   execute_on = timestep_end
  [../]
  [./slip_rate_10]
   type = MaterialStdVectorAux
   variable = slip_rate_10
   property = slip_increment
   index = 9
   execute_on = timestep_end
  [../]
  [./slip_rate_11]
   type = MaterialStdVectorAux
   variable = slip_rate_11   
   property = slip_increment
   index = 10
   execute_on = timestep_end
  [../]
  [./slip_rate_12]
   type = MaterialStdVectorAux
   variable = slip_rate_12
   property = slip_increment
   index = 11
   execute_on = timestep_end
  [../]
  [./slip_rate_13]
   type = MaterialStdVectorAux
   variable = slip_rate_13
   property = slip_increment
   index = 12
   execute_on = timestep_end
  [../]
  [./slip_rate_14]
   type = MaterialStdVectorAux
   variable = slip_rate_14
   property = slip_increment
   index = 13
   execute_on = timestep_end
  [../]
  [./slip_rate_15]
   type = MaterialStdVectorAux
   variable = slip_rate_15   
   property = slip_increment
   index = 14
   execute_on = timestep_end
  [../]
  [./slip_rate_16]
   type = MaterialStdVectorAux
   variable = slip_rate_16
   property = slip_increment
   index = 15
   execute_on = timestep_end
  [../]
  [./slip_rate_17]
   type = MaterialStdVectorAux
   variable = slip_rate_17
   property = slip_increment
   index = 16
   execute_on = timestep_end
  [../]
  [./slip_rate_18]
   type = MaterialStdVectorAux
   variable = slip_rate_18
   property = slip_increment
   index = 17
   execute_on = timestep_end
  [../]
  [./slip_rate_19]
   type = MaterialStdVectorAux
   variable = slip_rate_19   
   property = slip_increment
   index = 18
   execute_on = timestep_end
  [../]
  [./slip_rate_20]
   type = MaterialStdVectorAux
   variable = slip_rate_20
   property = slip_increment
   index = 19
   execute_on = timestep_end
  [../]
  [./slip_rate_21]
   type = MaterialStdVectorAux
   variable = slip_rate_21
   property = slip_increment
   index = 20
   execute_on = timestep_end
  [../]
  [./slip_rate_22]
   type = MaterialStdVectorAux
   variable = slip_rate_22
   property = slip_increment
   index = 21
   execute_on = timestep_end
  [../]
  [./slip_rate_23]
   type = MaterialStdVectorAux
   variable = slip_rate_23   
   property = slip_increment
   index = 22
   execute_on = timestep_end
  [../]
  [./slip_rate_24]
   type = MaterialStdVectorAux
   variable = slip_rate_24
   property = slip_increment
   index = 23
   execute_on = timestep_end
  [../]
  
  # build slip rate vector
  [./slip_rate_vector]
    type = BuildArrayVariableAux
    variable = slip_rate_vector
    component_variables = 'slip_rate_1 slip_rate_2 slip_rate_3 slip_rate_4 slip_rate_5 slip_rate_6 slip_rate_7 slip_rate_8 slip_rate_9 slip_rate_10 slip_rate_11 slip_rate_12 slip_rate_13 slip_rate_14 slip_rate_15 slip_rate_16 slip_rate_17 slip_rate_18 slip_rate_19 slip_rate_20 slip_rate_21 slip_rate_22 slip_rate_23 slip_rate_24'
  [../]
  
  [./edge_directional_derivative]
    type = ArrayDirectionalDerivative
    variable = dslip_increment_dedge
    gradient_variable = slip_rate_vector
    dislo_character = edge
  	execute_on = timestep_end
  [../]
  [./screw_directional_derivative]
    type = ArrayDirectionalDerivative
    variable = dslip_increment_dscrew
    gradient_variable = slip_rate_vector
    dislo_character = screw
  	execute_on = timestep_end
  [../]  
  
  [./elastic_energy]
    type = MaterialRealAux
    variable = elastic_energy
    property = elastic_energy
  [../]
  
  [./delastic_energy_dc]
    type = MaterialRealAux
    variable = delastic_energy_dc
    property = delastic_energy/dc
  [../]
  
  [./d2elastic_energy_dc2]
    type = MaterialRealAux
    variable = d2elastic_energy_dc2
    property = d^2elastic_energy/dc^2
  [../]
  
[]

[BCs]
  [./symmx]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0
  [../]
  [./symmy]
    type = DirichletBC
    variable = uy
    boundary = bottom
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
    variable = ux
    boundary = right
    function = '0.001*t*40'
  [../]
[]

[UserObjects]
  
  # poly_euler_angles.txt is a file with a list of Euler angles
  # on the structured mesh
  # row index: first growth along x, then y, then z
  # angles in degrees 
  [./prop_read]
    type = PropertyReadFile
    prop_file_name = 'EulerAngles.txt'
    nprop = 3
    read_type = block
    nblock = 8
    use_zero_based_block_indexing = false
  [../]
  
[]

[Transfers]
  [./from_c]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = fracture
    variable = c
    source_variable = c
  [../]
  
  [./to_elastic_energy]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = fracture
    variable = elastic_energy
    source_variable = elastic_energy
  [../]
  
  [./to_delastic_energy_dc]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = fracture
    variable = delastic_energy_dc
    source_variable = delastic_energy_dc
  [../]
  
  [./to_d2elastic_energy_dc2]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = fracture
    variable = d2elastic_energy_dc2
    source_variable = d2elastic_energy_dc2
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
    type = ComputeElasticityTensorCP
    C_ijkl = '2.3e5 1.35e5 1.35e5 2.3e5 1.35e5 2.3e5 1.17e5 1.17e5 1.17e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
  
  # coupled crystal plasticity and damage material
  # plastic_damage_prefactor from 0 to 1
  # set the fraction of plastic work that contributes to damage
  # the undamaged stress is used to calculate the plastic strain rate
  [./stress]
    type = CrystalPlasticityUndamagedStress
    crystal_plasticity_models = 'RPVsteel'
    tan_mod_type = exact
    maximum_substep_iteration = 1
    maxiter = 50000
    maxiter_state_variable = 50000
    rtol = 1e-3
    abs_tol = 1e-3
    c = c
    E_name = 'elastic_energy'
    D_name = 'degradation'
    use_current_history_variable = true
    use_snes_vi_solver = true
    plastic_damage_prefactor = 1.0
  [../]
  
  # model parameters calibrated for JRQ steel at room temperature
  [./RPVsteel]
    type = CrystalPlasticityFerriticSteel
    number_slip_systems = 24
    slip_sys_file_name = 'input_slip_sys_bcc24.txt'

	burgers_vector_mag = 0.000256
	shear_modulus = 86000.0
	a_self = 0.1
	a_col = 0.7
	
	ao = 0.042747355498052086
	xm = 0.14511893043165458
	
	# The ratio between Peierls' stress in different slip systems is kept as in:
	# Ghiath Monnet, Ludovic Vincent, Lionel Gelebart
    # Multiscale modeling of crystal plasticity in Reactor Pressure Vessel steels: Prediction of irradiation hardening
    # Journal of Nuclear Materials 514 (2019) 128-138
	const_slip_resistance_110 = 144.0
	const_slip_resistance_112_TW = 164.0
	const_slip_resistance_112_AT = 192.0
	
    k_0 = 0.5351213486387714
    y_c = 0.0016733349642994713
	init_rho_ssd = 2.6184166093059353
	
	init_rho_gnd_edge = 0.0
	init_rho_gnd_screw = 0.0
	rho_tol = 1.0
	dslip_increment_dedge = dslip_increment_dedge
	dslip_increment_dscrew = dslip_increment_dscrew
  [../]
  
  # degradation function for the positive part of the strain energy
  # it is used only in CrystalPlasticityUndamagedStress
  # therefore it must remain in the mechanics part
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '0.001'
    derivative_order = 2
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

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor -pc_hypre_boomeramg_print_statistics'
  petsc_options_value = 'hypre boomeramg 51 0.7 4 5 25 PMIS ext+i 2 0.3 0'

  line_search = 'none'
  
  nl_abs_tol = 1e-3
  nl_rel_tol = 1e-3

  nl_max_its = 10
  l_max_its = 50

  [./TimeStepper]
    type = ConstantDT
    dt = 0.01
    growth_factor = 10
  [../]

  dtmax = 0.01
  dtmin = 1.0e-40
  end_time = 0.01 # run until 100 s to see full crack propagation
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 1 #100
  [../]
  # uncomment for restart option
  #[./restart]
  #  type = Checkpoint
  #  num_files = 5
  #  interval = 10
  #[../]
[]
