# a variable growing from 0 to 1
# induces a linearly growing eigenstrain
# with nine independent components
# in this case non-zero components are xz, yz, zz

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz strain_xx strain_yy strain_zz strain_xy strain_xz strain_yz'
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

  # this variable changes from 0 to 1 with time
  # and induces a linearly increasing eigenstrain
  # with nine arbitrary components
  [./residual_def_level]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # AuxVariables can be introduced to visualize the components
  # of the residual deformation gradient
  [./residual_deformation_gradient_1_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'euler_angles.txt'
    nprop = 3
    read_type = element
  [../]
[]

[Functions]
  [./dts]
    type = PiecewiseConstant
    x = '0.0 1.0'
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

  [./fp_xx]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = plastic_deformation_gradient
	index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./fp_xy]
    type = RankTwoAux
    variable = fp_xy
    rank_two_tensor = plastic_deformation_gradient
	index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./fp_xz]
    type = RankTwoAux
    variable = fp_xz
    rank_two_tensor = plastic_deformation_gradient
	index_i = 0
    index_j = 2
    execute_on = timestep_end
  [../]
  [./fp_yx]
    type = RankTwoAux
    variable = fp_yx
    rank_two_tensor = plastic_deformation_gradient
	index_i = 1
    index_j = 0
    execute_on = timestep_end
  [../]
  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = plastic_deformation_gradient
	index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./fp_yz]
    type = RankTwoAux
    variable = fp_yz
    rank_two_tensor = plastic_deformation_gradient
	index_i = 1
    index_j = 2
    execute_on = timestep_end
  [../]
  [./fp_zx]
    type = RankTwoAux
    variable = fp_zx
    rank_two_tensor = plastic_deformation_gradient
	index_i = 2
    index_j = 0
    execute_on = timestep_end
  [../]
  [./fp_zy]
    type = RankTwoAux
    variable = fp_zy
    rank_two_tensor = plastic_deformation_gradient
	index_i = 2
    index_j = 1
    execute_on = timestep_end
  [../]
  [./fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = plastic_deformation_gradient
	index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]

  # variable increases from 0 to 1 at time 1
  [./residual_def_level]
    type = FunctionAux
    variable = residual_def_level
    function = 't'
    execute_on = timestep_begin
  [../]

  # AuxVariables can be introduced to visualize the components
  # of the residual deformation gradient
  [./residual_deformation_gradient_1_zz]
    type = RankTwoAux
    variable = residual_deformation_gradient_1_zz
    rank_two_tensor = residual_deformation_gradient_1
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]

[]

# the element is constrained in all directions on the back surface
# and eigenstrains are visible as displacement on the front surface
# there is no displacement applied, the only displacement
# visible is due to the eigenstrain
[BCs]
  [./z0_back]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./y0_back]
    type = DirichletBC
    variable = disp_y
    boundary = back
    value = 0.0
  [../]

  [./x0_back]
    type = DirichletBC
    variable = disp_x
    boundary = back
    value = 0.0
  [../]
[]

# parameters for additively manufactured stainless steel
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
    eigenstrain_names = 'residual_eigenstrain_1'
    tan_mod_type = exact
    maximum_substep_iteration = 2
	  maxiter = 500
	  maxiter_state_variable = 500
  [../]
  # this eigenstrain is a method to introduce initial residual stress
  # the components correspond to eigenstrain components in the following order
  # 00, 10, 20, 01, 11, 21, 02, 12, 22
  # in this case only xz, yz, zz are non-zero
  # this are components of a residual deformation gradient
  # not of a symmetric strain tensor
  # residual_deformation_gradient_1 is useful for visualization
  # this is just one element without external displacement BC
  # therefore the stress is zero even though there is a non-zero displacement
  # by definition of eigenstrain
  [residual_eigenstrain_1]
    type = ComputeCrystalPlasticityResidualEigenstrain
    eigenstrain_name = residual_eigenstrain_1
    deformation_gradient_name = residual_deformation_gradient_1
    residual_def_level = residual_def_level
    residual_def_components = '0.0 0.0 0.0 0.0 0.0 0.0 0.01 0.01 0.01'
  []
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
	  y_c = 0.0019545318633428007
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

  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8

  [./TimeStepper]
    type = FunctionDTGrowth
    function = dts
	  cutback_factor_at_failure = 0.5
	  growth_factor = 1.1
  [../]

  end_time = 1.0
  dtmin = 0.000001

[]

[Outputs]
  exodus = true
  csv = false
  interval = 5
[]
