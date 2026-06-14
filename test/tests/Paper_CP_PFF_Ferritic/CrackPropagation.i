# BCC model for irradiated RPV steel
# JRQ ferritic steel model parameters
# crystal plasticity with phase field fracture
# pre-existing crack to show damage propagation
# It corresponds to the simulation in Figure 8 of
# Michael Salvini, Nicolò Grilli, Mahmoud Mostafavi, Guilherme Correa Soares, Matti Lindroos, Marta Serrano, Christopher Truman
# Understanding the Effect of Neutron Dose on Plastic Deformation, Fracture Nucleation, 
# and Toughness Using a Crystal Plasticity Phase Field Fracture Model
# and a Microscale J-Integral methodology
# https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6834348

[GlobalParams]
  displacements = 'ux uy uz'
[]

[Mesh]
  [./generated_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    ymin = 0.0
    xmax = 250.305
    xmin = 0.0
    ymax = 201.583
    zmin = 0.0
    zmax = 2.0
    nx = 128
    ny = 100
    nz = 1
    elem_type = HEX8
  [../]

  # boundary to apply displacement control
  # on the left hand side of the crack
  [./bottom_left_edge]
    type = BoundingBoxNodeSetGenerator
    new_boundary = 'bl_crack'
    input = 'generated_mesh'
    top_right = '125.6 0.5 3'
    bottom_left = '-1 -0.5 -1'
  [../]
[]

[Functions]
  [./pull_right_x]
    type = ParsedFunction
    expression = '(t / 4)'
  [../]
[]

[ICs]
  # define pre-existing crack using the fracture phase field
  [./ic_c]
    type = FunctionIC
    variable = c
    function = '1.0 * if( y >= 40.0 , exp(-abs(((y-40.0)^2) + ((x-125.575)^2)) / 50.0) , 1.0) * exp(-abs((x-125.575)^2) / 20.0)'
  [../]
[]

[AuxVariables]

  [./elastic_energy]
   order = FIRST
   family = MONOMIAL
  [../]

  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
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

[]

[Physics/SolidMechanics/QuasiStatic/all]
  strain = FINITE
  add_variables = true
  generate_output = 'stress_xx stress_xy stress_yy stress_yz stress_zz stress_xz strain_xx strain_xy strain_yy strain_yz strain_zz strain_xz'
[]

[Modules/PhaseField/Nonconserved/c]
  free_energy = F
  kappa = kappa_op
  mobility = L
[]

# off diagonal Jacobian kernels
[Kernels]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = ux
    component = 0
    c = c
    use_displaced_mesh = true
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = uy
    component = 1
    c = c
    use_displaced_mesh = true
  [../]
  [./solid_z]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = uz
    component = 2
    c = c
    use_displaced_mesh = true
  [../]
  
  [./off_disp]
    type = AllenCahnElasticEnergyOffDiag
    variable = c
    displacements = 'ux uy uz'
    mob_name = L
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]

  [./elastic_energy]
    type = MaterialRealAux
    variable = elastic_energy
    property = elastic_energy
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
[]

# bounds used to impose that damage can only grow
[Bounds]
  [./irreversibility]
    type = VariableOldValueBounds
    variable = bounds_dummy
    bounded_variable = c
    bound_type = lower
  [../]
  [./upper]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = c
    bound_type = upper
    bound_value = 1.0
  [../]
[]

[BCs]

  [./edge_damage_l]
    type = DirichletBC
    variable = c
    boundary = left
    value = 0
  [../]

  [./edge_damage_r]
    type = DirichletBC
    variable = c
    boundary = right
    value = 0
  [../]

  [./bl_crack_no_y_move]
    type = DirichletBC
    variable = uy
    boundary = bl_crack
    value = 0
  [../]

  [./left]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0
  [../]

  [./back]
    type = DirichletBC
    variable = uz
    boundary = back
    value = 0
  [../]

  [./right_x]
    type = FunctionDirichletBC
    variable = ux
    boundary = right
    function = pull_right_x
  [../]
[]

[UserObjects]
  
  # poly_euler_angles.txt is a file with a list of Euler angles
  # on the structured mesh
  # row index: first growth along x, then y, then z
  # angles in degrees 
  [./prop_read]
    type = PropertyReadFile
    prop_file_name = 'EBSD.txt'
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
    type = ComputeElasticityTensorCP
    C_ijkl = '2.3e5 1.35e5 1.35e5 2.3e5 1.35e5 2.3e5 1.17e5 1.17e5 1.17e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]

  [./stress]
    type = CrystalPlasticityUndamagedStress
    crystal_plasticity_models = 'irradiatedRPVsteel'
    tan_mod_type = exact
    maximum_substep_iteration = 1
	  rtol = 1e-5
	  abs_tol = 1e-5
	  c = c
    E_name = 'elastic_energy'
    D_name = 'degradation'
    use_current_history_variable = true
    plastic_damage_prefactor = 1.0 # plastic work contribution to damage
  [../]
  
  # model parameters calibrated for JRQ ferritic steel at room temperature
  [./irradiatedRPVsteel]
    type = CrystalPlasticityFerriticSteel
    number_slip_systems = 24
    slip_sys_file_name = 'input_slip_sys_bcc24.txt'

#   uncomment to activate neutron irradiation at 1.04 dpa
#   see table 4 in the paper 
#   is_irradiated = true
#   annihilate_DL_on_all_slip_systems = true
#   init_C_DL = 624.0
#   C_DL_diameter = 0.0025
#   a_DL = 1.0
#   lambda_DL = 1.0
#   init_C_SC = 59300.0
#   C_SC_diameter = 0.00154
#   a_SC = 0.10
#   lambda_SC = 1.0

#   see table 3 in the paper
    burgers_vector_mag = 0.000256
    shear_modulus = 86000.0
    a_self = 0.1
    a_col = 0.7
    ao = 0.001
    xm = 0.05

#   see table 4 in the paper
    const_slip_resistance_110 = 90.0
    const_slip_resistance_112_TW = 105.0
    const_slip_resistance_112_AT = 120.0
    k_0 = 0.56
    y_c = 0.00398
    init_rho_ssd = 4.32

    init_rho_gnd_edge = 0.0
    init_rho_gnd_screw = 0.0
    rho_tol = 1.0
    dslip_increment_dedge = dslip_increment_dedge
    dslip_increment_dscrew = dslip_increment_dscrew
  [../]

  [./degradation]
    type = DerivativeParsedMaterial
    property_name = degradation
    coupled_variables = 'c'
    expression = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '1.0e-2'
    derivative_order = 2
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '10 4.0 1e-3'
  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    property_name = L
    expression = '1.0/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    property_name = kappa_op
    expression = 'gc_prop * l'
  [../]
  [./local_fracture_energy]
    type = DerivativeParsedMaterial
    property_name = local_fracture_energy
    coupled_variables = 'c'
    material_property_names = 'gc_prop l'
    expression = 'c^2 * gc_prop / 2 / l'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    coupled_variables = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    property_name = F
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

  petsc_options_iname = '-pc_type  -snes_type'
  petsc_options_value = 'lu vinewtonrsls'
  
  nl_abs_tol = 1e-5
  nl_rel_tol = 1e-5
  
  dtmax = 0.005
  dtmin = 0.000000001

  nl_max_its = 50
  l_max_its = 100

  start_time = 0.0
  end_time = 500.0

  [./TimeStepper]
    type = ConstantDT
    dt = 0.005
    growth_factor = 2.0
  [../]
[]

[Outputs]
  [./out]
    type = Exodus
    time_step_interval = 1000
  [../]
[]
