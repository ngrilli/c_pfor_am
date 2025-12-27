# damage growth induced by positive volumetric strain
# crack propagation in a 3D geometry with pre-existing crack
# modelled as an initial phase field

# Free energy is decomposed into volumetric and non-volumetric parts as in:
# Nicolo Grilli and Marisol Koslowski
# The effect of crystal anisotropy and plastic response
# on the dynamic fracture of energetic materials
# Journal of Applied Physics 126, 155101 (2019)

# Bicrystal: two crystals on left and right
# cleavage plane for fracture determined by slip plane normal
# of first slip system (slip_sys_index = 0)
# therefore crack is expected to deviate at the
# grain boundary

# Bicrystal on the right is 90 degrees rotated compared
# with the one on the left, therefore crack deviates
# by about 90 degrees on the right hand side

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [dual_block]
    type = CartesianMeshGenerator
    dim = 3
    dx = '5 5'
    dy = '10.0'
    dz = '0.5'
    ix = '10 10'
    iy = '20'
    iz = '1'
    subdomain_id = '1 2'
  []
[]

[AuxVariables]
  [./pk2_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pk2_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pk2_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hist]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./bounds_dummy]
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
  [./tdisp]
    type = ParsedFunction
    value = '0.05*t'
  [../]
[]

# c variable must be added manually
# because phase field action is not used 
[Variables]
  [./c]
    family = LAGRANGE
    order = FIRST
  [../]
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
[]

[ICs]
  [./ic_c]
    type = FunctionIC
    variable = c
    function = '0.98*if((x-0.1)/5.0,0.0,1.0)*exp(-abs(x+y-9.0))'
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

[Kernels]
  # kernels for phase field damage
  [./ACbulk]
    type = AllenCahn
    variable = c
    f_name = F
  [../]
  [./ACInterfaceCleavageFracture]
    type = ACInterfaceSlipPlaneFracture
    variable = c
    beta_penalty = 10
    slip_sys_index = 0 # slip system index whose plane normal is the cleavage plane
  [../]
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]
  
  # off diagonal Jacobian kernels
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
    use_displaced_mesh = true
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
    use_displaced_mesh = true
  [../]
  [./solid_z]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_z
    component = 2
    c = c
    use_displaced_mesh = true
  [../]
  [./off_disp]
    type = AllenCahnElasticEnergyOffDiag
    variable = c
    displacements = 'disp_x disp_y disp_z'
    mob_name = L
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./pk2_zz]
   type = RankTwoAux
   variable = pk2_zz
   rank_two_tensor = pk2
   index_j = 2
   index_i = 2
   execute_on = timestep_end
  [../]
  [./pk2_yy]
   type = RankTwoAux
   variable = pk2_yy
   rank_two_tensor = pk2
   index_j = 1
   index_i = 1
   execute_on = timestep_end
  [../]
  [./pk2_xx]
   type = RankTwoAux
   variable = pk2_xx
   rank_two_tensor = pk2
   index_j = 0
   index_i = 0
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
  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = fp
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
  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = lage
    index_j = 2
    index_i = 2
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
  [./e_xx]
    type = RankTwoAux
    variable = e_xx
    rank_two_tensor = lage
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./hist]
    type = MaterialRealAux
    variable = hist
    property = hist
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
  [./disp_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./disp_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./disp_back]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  [../]
  [./disp_top]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = 'tdisp'
  [../]
[]

[UserObjects]
  [./slip_rate_gss_1]
    type = CrystalPlasticitySlipRateCleavage
    variable_size = 12
    slip_sys_file_name = 'input_slip_sys.txt'
    num_slip_sys_flowrate_props = 2
    flowprops = '1 12 0.001 0.14'
    uo_state_var_name = state_var_gss_1
    block = '1'
  [../]
  [./slip_rate_gss_2]
    type = CrystalPlasticitySlipRateCleavage
    variable_size = 12
    slip_sys_file_name = 'input_slip_sys.txt'
    num_slip_sys_flowrate_props = 2
    flowprops = '1 12 0.001 0.14'
    uo_state_var_name = state_var_gss_2
    block = '2'
  [../]
  [./slip_resistance_gss_1]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 12
    uo_state_var_name = state_var_gss_1
    block = '1'
  [../]
  [./slip_resistance_gss_2]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 12
    uo_state_var_name = state_var_gss_2
    block = '2'
  [../]
  [./state_var_gss_1]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 3 6 12'
    group_values = '420 420 420'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss_1
    scale_factor = 1.0
    block = '1'
  [../]
  [./state_var_gss_2]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 3 6 12'
    group_values = '420 420 420'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss_2
    scale_factor = 1.0
    block = '2'
  [../]
  [./state_var_evol_rate_comp_gss_1]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 12
    hprops = '1 631.2 462 2.5'
    uo_slip_rate_name = slip_rate_gss_1
    uo_state_var_name = state_var_gss_1
    block = '1'
  [../]
  [./state_var_evol_rate_comp_gss_2]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 12
    hprops = '1 631.2 462 2.5'
    uo_slip_rate_name = slip_rate_gss_2
    uo_state_var_name = state_var_gss_2
    block = '2'
  [../]
[]

[Materials]
  [./crysp_1] # class outputting slip plane for the cleavage plane 
    type = FiniteStrainUObasedCPCleavage
    c = c
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss_1'
    uo_slip_resistances = 'slip_resistance_gss_1'
    uo_state_vars = 'state_var_gss_1'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss_1'
    E_name = 'elastic_energy'
    D_name = 'degradation'
    use_current_history_variable = false
    maxiter = 1000
    maximum_substep_iteration = 5
    use_line_search = true
    block = '1'
  [../]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.624e5 .92e5 .69e5 1.624e5 .69e5 1.807e5 0.467e5 0.467e5 0.352e5'
    fill_method = symmetric9
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    block = '1'
  [../]
  [./crysp_2] # new class with damage (volumetric)
    type = FiniteStrainUObasedCPCleavage
    c = c
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss_2'
    uo_slip_resistances = 'slip_resistance_gss_2'
    uo_state_vars = 'state_var_gss_2'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss_2'
    E_name = 'elastic_energy'
    D_name = 'degradation'
    use_current_history_variable = false
    maxiter = 1000
    maximum_substep_iteration = 5
    use_line_search = true
    block = '2'
  [../]
  [./elasticity_tensor_2]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.624e5 .92e5 .69e5 1.624e5 .69e5 1.807e5 0.467e5 0.467e5 0.352e5'
    fill_method = symmetric9
    euler_angle_1 = 90.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    block = '2'
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    property_name = degradation
    coupled_variables = 'c'
    expression = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '5.0e-4'
    derivative_order = 2
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '1.0 1.0 4e-3'
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
  dt = 0.005
  solve_type = 'PJFNK'

  # variational inequalities option vinewtonrsls
  # is used to impose that damage can only grow
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'

  # high tolerance to make the test faster
  # keep default tolerance for production simulations  
  nl_rel_tol = 1e-3 #1e-8
  
  l_max_its = 20
  nl_max_its = 20

  line_search='none'

  dtmax = 10.0
  dtmin = 0.0000001

  num_steps = 1 #1000 to see full crack propagation
[]

[Outputs]
  exodus = true
  interval = 1 #5
[]
