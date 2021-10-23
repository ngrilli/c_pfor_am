# damage growth induced by positive volumetric strain
# Free energy is decomposed into volumetric and non-volumetric parts as in:
# Nicolo Grilli and Marisol Koslowski
# The effect of crystal anisotropy and plastic response 
# on the dynamic fracture of energetic materials
# Journal of Applied Physics 126, 155101 (2019)

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
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
  [./fp_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_increment]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
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
  [./off_disp]
    type = AllenCahnElasticEnergyOffDiag
    variable = c
    displacements = 'disp_x disp_y'
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
  [./fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = fp
    index_j = 2
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
  [./slip_inc]
   type = MaterialStdVectorAux
   variable = slip_increment
   property = slip_rate_gss
   index = 0
   execute_on = timestep_end
  [../]
  [./gss]
    type = MaterialStdVectorAux
    variable = gss
    property = state_var_gss
    index = 0
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
  [./disp_right]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0
  [../]
  [./disp_back]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  [../]
  [./disp_front]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = '0.01*t'
  [../]
[]

[UserObjects]
  [./slip_rate_gss]
    type = CrystalPlasticitySlipRateGSS
    variable_size = 12
    slip_sys_file_name = input_slip_sys.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1'
    uo_state_var_name = state_var_gss
  [../]
  [./slip_resistance_gss]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 12
    uo_state_var_name = state_var_gss
  [../]
  [./state_var_gss]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 4 8 12'
    group_values = '60.8 60.8 60.8'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 12
    hprops = '1.0 541.5 109.8 2.5'
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  [../]
[]

[Materials]
  [./crysp] # new class with damage (volumetric)
    type = FiniteStrainUObasedCPDamageVol
	c = c
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
    E_name = 'elastic_energy'
    D_name = 'degradation'
    use_current_history_variable = true
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '1.0e-6'
    derivative_order = 2
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '2.0 0.5 1e-6'
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
  dt = 0.001
  solve_type = 'PJFNK'

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomerang
  line_search='none'
  nl_abs_tol = 1e-8
  nl_rel_step_tol = 1e-8
  dtmax = 10.0
  nl_rel_tol = 1e-8
  dtmin = 0.001
  
  # select 300 steps to see elasto-plastic curve
  # followed by damage
  num_steps = 2 #300
  nl_abs_step_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
