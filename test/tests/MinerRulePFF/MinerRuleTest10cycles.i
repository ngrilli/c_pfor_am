# Miner's rule one element test
# 10 load cycles representing 10000 cycles

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmax = 1
  ymax = 1
  zmax = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules]
  [./PhaseField]
    [./Nonconserved]
      [./c]
        free_energy = F
        kappa = kappa_op
        mobility = L
      [../]
    [../]
  [../]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz strain_xx strain_yy strain_zz strain_xy strain_xz strain_yz'
      [../]
    [../]
  [../]
[]

[ICs]
[]

[Functions]
  # one cycle per second
  # strain amplitude corresponding to a 100 MPa stress
  [./pull]
    type = ParsedFunction
    value = '0.001638 * sin(6.28*t)'
  [../]
[]

[AuxVariables]
  [./alpha_cyclic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fatigue_degradation]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  [./NS_curve]
    order = CONSTANT
    family = MONOMIAL  
  [../]
  [./cyclic_stress_history]
    order = CONSTANT
    family = MONOMIAL 
  [../]
[]

[Kernels]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
  [../]
  [./solid_z]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_z
    component = 2
    c = c
  [../]
  [./off_disp]
    type = AllenCahnElasticEnergyOffDiag
    variable = c
    displacements = 'disp_x disp_y disp_z'
    mob_name = L
  [../]
[]

[AuxKernels]
  [./alpha_cyclic]
    type = MaterialRealAux
    variable = alpha_cyclic
    property = alpha_cyclic 
  [../]
  [./fatigue_degradation]
    type = MaterialRealAux
    variable = fatigue_degradation
    property = fatigue_degradation
  [../]  
  [./NS_curve]
    type = MaterialRealAux
    variable = NS_curve
    property = NS_curve
  [../]
  [./cyclic_stress_history]
    type = MaterialRealAux
    variable = cyclic_stress_history
    property = cyclic_stress_history
  [../]
[]

[BCs]
  [./block_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./block_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./block_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  [./x_load]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = pull
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '30000000 4 1e-4'
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
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.26
    youngs_modulus = 61040000000 # Pa
  [../]
  
  # tau_cyclic_stress_history much larger than the load cycle period
  # so that cyclic_stress_history is not underestimated
  # and fatigue life overpredicted
  [./damage_stress]
    type = ComputeLinearElasticPFFractureCyclic
    c = c
    E_name = 'elastic_energy'
    D_name = 'degradation'
    F_name = 'local_fracture_energy'
    decomposition_type = strain_spectral
    cycles_per_unit_time = 1000 # cycles per second
    tau_cyclic_stress_history = 100.0 # second
    fatigue_degradation_name = fatigue_degradation
    NS_curve_name = NS_curve
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    property_name = degradation
    coupled_variables = 'c'
    expression = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '0.001'
    derivative_order = 2
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
    coupled_variables = 'c'
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    property_name = F
  [../]
  
  # Number of cycles to failure at a given stress level
  # calculated based on cyclic_stress_history
  [./NS_curve]
    type = DerivativeParsedMaterial
    material_property_names = 'cyclic_stress_history'
    property_name = NS_curve
    expression = '10^(8 - 0.04e-6 * cyclic_stress_history)' # N as a function of S, thus invert the standard SN plot
  [../]
  
  # this is effectively a reduction of gc_prop
  # in this case, cyclic degradation starts when the 
  # Miner rule variable goes above 0.9
  [./fatigue_degradation]
    type = DerivativeParsedMaterial
    material_property_names = 'alpha_cyclic'
    property_name = fatigue_degradation
    expression = 'if(alpha_cyclic/1.8,10*(1-alpha_cyclic),1)'
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

  solve_type = PJFNK
  
  petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor'
  petsc_options_value = 'hypre boomeramg 0.7 4 5 25 HMIS ext+i 2 0.3'

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-5
  
  l_max_its = 100
  nl_max_its = 20

  dt = 0.004
  dtmin = 1e-9
  end_time = 0.004 # run until time 14 to see cyclic fatigue failure after about 10 cycles
[]

[Outputs]
  exodus = true
  interval = 5
[]
