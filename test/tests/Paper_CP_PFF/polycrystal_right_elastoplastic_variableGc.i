# This simulation corresponds to figure 11b in:
# M Salvini, N Grilli, E Demir, S He, T Martin, P Flewitt, M Mostafavi, C Truman, D Knowles
# Effect of grain boundary misorientation and carbide precipitation on damage initiation: 
# A coupled crystal plasticity and phase field damage study
# International Journal of Plasticity 172 (2024) 103854
# https://www.sciencedirect.com/science/article/pii/S0749641923003388

[GlobalParams]
  displacements = 'ux uy uz'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = 0.0
  xmax = 1020.0
  ymin = 0.0
  ymax = 1014.0
  zmin = 0.0
  zmax = 5.0
  nx = 170
  ny = 169
  nz = 1
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

[Functions]

  [pull_y]
    type = PiecewiseLinear
    data_file = 'Right_Slip_Creep.csv'
    format = columns
  []

  # Plastic deformation followed by creep by changing the slip rate prefactor with time  
  [./creep_rate_prefactor]
    type = PiecewiseConstant
    x = '0.0 1.4'
    y = '0.001 3.0e-8'
  [../]

[]

[Modules/TensorMechanics/Master/all]
  add_variables = true
  strain = FINITE
  generate_output = 'strain_xx strain_xy strain_xz strain_yy strain_yz strain_zz stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz'
[]

[Modules/PhaseField/Nonconserved/c]
  free_energy = F
  kappa = kappa_op
  mobility = L
[]

# initial condition for the damage
[ICs]
  [./ic_c]
    type = FunctionIC
    variable = c
    function = '0.0'
  [../]
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

  [./plastic_work]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./bounds_dummy]
  [../]
  
  # free energy terms that contribute to damage
  [./elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  # Gc visualization
  [./Gc]
    order = CONSTANT
    family = MONOMIAL
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

  [./plastic_work]
    type = MaterialRealAux
    variable = plastic_work
    property = plastic_work
    execute_on = timestep_end
  [../]
  
  [./elastic_energy]
    type = MaterialRealAux
    variable = elastic_energy
    property = elastic_energy
    execute_on = timestep_end
  [../]
  
  [./Gc]
    type = MaterialRealAux
    variable = Gc
    property = gc_prop
    execute_on = timestep_end
  [../]

[]

# bounds used to impose that damage can only grow
[Bounds]
  [./irreversibility]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = lower
  [../]
  [./upper]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = upper
    bound_value = 1.0
  [../]
[]

# tension boundary conditions
# lateral boundaries blocked
[BCs]
  [y_pull_top]
    type = FunctionDirichletBC
    variable = uy
    boundary = top
    function = pull_y
  []
  [x_left]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0.0
  []
  [y_bottom]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0.0
  []
  [z_back]
    type = DirichletBC
    variable = uz
    boundary = back
    value = 0.0
  []
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPGrain
    C_ijkl = '2.046e5 1.377e5 1.377e5 2.046e5 1.377e5 2.046e5 1.262e5 1.262e5 1.262e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
  
  # plastic_damage_prefactor from 0 to 1
  # set the fraction of plastic work that contributes to damage
  [./stress]
    type = ComputeCrystalPlasticityStressDamage
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    maximum_substep_iteration = 2
	maxiter = 500
	maxiter_state_variable = 500
	c = c
    E_name = 'elastic_energy'
    D_name = 'degradation'
    use_current_history_variable = true
    plastic_damage_prefactor = 1.0
  [../]
  
  # 316H stainless steel parameters
  [./trial_xtalpl]
    type = CrystalPlasticityDislocationUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
	ao = 0.0 # creep term acts as slip before time = 1.4
	xm = 0.1
	creep_xm = 0.1
	creep_ao_function = creep_rate_prefactor
	burgers_vector_mag = 0.000256
	shear_modulus = 86000 # MPa
	alpha_0 = 0.3
	r = 1.4
	tau_c_0 = 0.112
	k_0 = 1.13
	y_c = 0.0013
	init_rho_ssd = 4.921
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
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '1.0e-3'
    derivative_order = 2
  [../]
  
  # gc_prop is the fracture energy
  # if plastic work is included
  # fracture takes place at lower strain
  # l is the characteristic crack diffusion length
  # ideally element size should be about one fourth of l
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'l visco'
    prop_values = '3.0 1e-3'
  [../]
  
  # gc_prop is read from file
  # Gc is a function of coordinates
  [./Gc_from_file_object]
    type = FileParsedMaterial
    prop_name = 'gc_prop'
    read_prop_user_object = Gc_read
  [../]
  
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l'
  [../]
  [./local_fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'c'
    material_property_names = 'gc_prop l'
    function = 'c^2 * gc_prop / 2 / l'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    args = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    f_name = F
  [../]
[]

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'euler_angles.txt'
    nprop = 3
    read_type = element
  [../]
  
  [./Gc_read]
    type = PropertyReadFile
    prop_file_name = 'GB_Gc_078_linear.txt'
    nprop = 1
    read_type = element
  [../]
[]

[Postprocessors]
[]

[Executioner]
  type = Transient
  
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type  -snes_type'
  petsc_options_value = 'lu vinewtonrsls'

  line_search = 'none'

  automatic_scaling = true
  
  nl_max_its = 20
  l_max_its = 50
  
  nl_rel_tol = 1e-8
  
  [TimeStepper]
    type = IterationAdaptiveDTMax
    optimal_iterations = 25
    iteration_window = 15
    growth_factor = 1.25
  	reject_large_step = True
    cutback_factor = 0.5
    dt = 0.0005
	upper_limit_dt = 0.0005
  []
  
  end_time = 0.005 # run until 12503.8
  dtmin = 1.0e-20
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Outputs]
  [./out]
    type = Exodus
    sync_times = '0.005 0.1 0.5 1.0 1.4 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0 10.5 11.0 11.5 12.0 12.5 13.0 13.5 14.0 14.5 15.0 50 150 500 1000 1500 12503.8'
	sync_only = true
  [../]
[]
