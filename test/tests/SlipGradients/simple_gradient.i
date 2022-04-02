# simple test with only auxiliary variables 
# defined by function that are linear in space
# the Euler angles are such that the slip direction
# of the first slip system is along the x axis
# while the slip plane normal is along the z axis
# therefore you can change the order of the auxvariables in
# build_var_vector and observe that the directional derivatives
# of the first auxvariable are taken exactly along the x and y direction
# depending if you look at d_vettore_edge or d_vettore_screw

[GlobalParams]
  displacements = 'ux uy uz'
[]

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
  displacements = 'ux uy uz'
[]

[Variables]
  [./ux]
    order = FIRST
    family = LAGRANGE
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
  generate_output = stress_xz
[]

[AuxVariables]

  [./variabile_ausiliaria_x]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./variabile_ausiliaria_y]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./variabile_ausiliaria_z]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./zero_first]
    order = FIRST
    family = MONOMIAL
  [../]

  [./vettore]
    order = FIRST
    family = MONOMIAL
    components = 12
  [../]

  [./d_vettore_edge]
    order = CONSTANT
    family = MONOMIAL
    components = 12
  [../]
  
  [./d_vettore_screw]
    order = CONSTANT
    family = MONOMIAL
    components = 12
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

[Functions]

  [./funzione_x]
    type = ParsedFunction
    value = 'x'
  [../]
  
  [./funzione_y]
    type = ParsedFunction
    value = 'y'
  [../]
  
  [./funzione_z]
    type = ParsedFunction
    value = 'z'
  [../]

[]

[AuxKernels]

  [./variabile_ausiliaria_x]
    type = FunctionAux
    variable = variabile_ausiliaria_x   
    function = funzione_x
    execute_on = timestep_begin
  [../]
  
  [./variabile_ausiliaria_y]
    type = FunctionAux
    variable = variabile_ausiliaria_y
    function = funzione_y
    execute_on = timestep_begin
  [../]
  
  [./variabile_ausiliaria_z]
    type = FunctionAux
    variable = variabile_ausiliaria_z 
    function = funzione_z
    execute_on = timestep_begin
  [../]

  [./edge_directional_derivative]
    type = ArrayDirectionalDerivative
    variable = d_vettore_edge
    gradient_variable = vettore
    dislo_character = edge
  	execute_on = timestep_end
  [../]
  
  [./screw_directional_derivative]
    type = ArrayDirectionalDerivative
    variable = d_vettore_screw
    gradient_variable = vettore
    dislo_character = screw
  	execute_on = timestep_end
  [../]

  # move the non-zero variables at different positions
  # of the vector to see the gradient at different positions
  # in d_vettore_edge and d_vettore_screw
  # each auxvariable will be derived with respect to
  # corresponding slip system directions
  # depending on the position in the vector 
  [./build_var_vector]
    type = BuildArrayVariableAux
    variable = vettore
    component_variables = 'variabile_ausiliaria_y variabile_ausiliaria_x zero_first zero_first zero_first variabile_ausiliaria_z zero_first zero_first zero_first zero_first zero_first zero_first'
  	execute_on = timestep_end
  [../]
  
[]

[BCs]
  [./back_y]
    type = DirichletBC
    variable = uy
    boundary = back
    value = 0
  [../]
  [./back_x]
    type = DirichletBC
    variable = ux
    boundary = back
    value = 0
  [../]
  [./back_z]
    type = DirichletBC
    variable = uz
    boundary = back
    value = 0
  [../]
  [./shear_load_x]
    type = FunctionDirichletBC
    variable = ux
    boundary = front
    function = 't'
  [../]
  [./front_y]
    type = DirichletBC
    variable = uy
    boundary = front
    value = 0  
  [../]
  [./front_z]
    type = DirichletBC
    variable = uz
    boundary = front
    value = 0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPGrain
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
  [./stress]
    type = ComputeDislocationCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    maximum_substep_iteration = 1
  [../]
  [./trial_xtalpl]
    type = CrystalPlasticityDislocationUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
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
  dt = 0.001
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
  nl_abs_tol = 1e-8
  nl_rel_step_tol = 1e-8
  dtmax = 0.1
  nl_rel_tol = 1e-8
  end_time = 1
  dtmin = 0.00001
  num_steps = 1 #100 to see plastic shear stress
  nl_abs_step_tol = 1e-8
[]

[Outputs]
  exodus = true
[]
