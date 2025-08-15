# example with two sensitivities

vol_frac = 0.4
E0 = 1e5
Emin = 1e-4
power = 2

[GlobalParams]
[]

[Mesh]
  [MeshGenerator]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 24
    ny = 12
    nz = 12
    xmin = 0
    xmax = 20
    ymin = 0
    ymax = 10
    zmin = 0
    zmax = 10
  []
  [middle_bottom_left_edge]
    type = ExtraNodesetGenerator
    input = MeshGenerator
    new_boundary = pull
    coord = '0 5 5'
  []
[]

[Variables]
  [Dc]
    initial_condition = -1.0
  []

  # displacement for elastic problem  
  [./ue_x]
  [../]
  [./ue_y]
  [../]
  [./ue_z]
  [../]
  
  # displacement for thermo-elastic problem  
  [./ut_x]
  [../]
  [./ut_y]
  [../]
  [./ut_z]
  [../]
[]

[AuxVariables]
  [sensitivity]
    family = MONOMIAL
    order = FIRST
    initial_condition = -1.0
    [AuxKernel]
      type = MaterialRealAux
      variable = sensitivity
      property = mech_sensitivity
      execute_on = LINEAR
    []
  []

  [compliance]
    family = MONOMIAL
    order = CONSTANT
  []

  [mat_den]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${vol_frac}
  []
  
  [Dc_elem]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = -1.0
    [AuxKernel]
      type = SelfAux
      variable = Dc_elem
      v = Dc
      execute_on = 'TIMESTEP_END'
    []
  []
[]

[Kernels]
  [diffusion]
    type = FunctionDiffusion
    variable = Dc
    function = 0.05
  []
  [potential]
    type = Reaction
    variable = Dc
  []
  [source]
    type = CoupledForce
    variable = Dc
    v = sensitivity
  []
  
  # for mechanics problem, use kernels instead of action
  # so that two problems with different displacement variables
  # can be solved at the same time, for example these are the 3 equations
  # for equilibrium for the mechanics problem
  [./ue_x]
    type = StressDivergenceTensors
    variable = ue_x
    component = 0
    base_name = "mech"
    displacements = 'ue_x ue_y ue_z'
  [../]
  [./ue_y]
    type = StressDivergenceTensors
    variable = ue_y
    component = 1
    base_name = "mech"
    displacements = 'ue_x ue_y ue_z'
  [../]
  [./ue_z]
    type = StressDivergenceTensors
    variable = ue_z
    component = 2
    base_name = "mech"
    displacements = 'ue_x ue_y ue_z'
  [../]
  
  # and these other 3 equations for equilibrium for the thermo-mechanical problem
  [./ut_x]
    type = StressDivergenceTensors
    variable = ut_x
    component = 0
    base_name = "mech"
    displacements = 'ut_x ut_y ut_z'
    eigenstrain_names = eigenstrain
  [../]
  [./ut_y]
    type = StressDivergenceTensors
    variable = ut_y
    component = 1
    base_name = "mech"
    displacements = 'ut_x ut_y ut_z'
    eigenstrain_names = eigenstrain
  [../]
  [./ut_z]
    type = StressDivergenceTensors
    variable = ut_z
    component = 2
    base_name = "mech"
    displacements = 'ut_x ut_y ut_z'
    eigenstrain_names = eigenstrain
  [../]
[]

[BCs]
  [no_e_x]
    type = DirichletBC
    variable = ue_x
    boundary = right
    value = 0.0
  []
  [no_e_y]
    type = DirichletBC
    variable = ue_y
    boundary = right
    value = 0.0
  []
  [no_e_z]
    type = DirichletBC
    variable = ue_z
    boundary = right
    value = 0.0
  []
  
  [boundary_penalty]
    type = ADRobinBC
    variable = Dc
    boundary = 'left top front back'
    coefficient = 10
  []
  
  # boundary condition for substrate in thermal simulation
  [no_t_x]
    type = DirichletBC
    variable = ut_x
    boundary = 'right left'
    value = 0.0
  []
  [no_t_y]
    type = DirichletBC
    variable = ut_y
    boundary = 'right left'
    value = 0.0
  []
  [no_t_z]
    type = DirichletBC
    variable = ut_z
    boundary = 'right left'
    value = 0.0
  []
[]

[NodalKernels]
  [pull]
    type = NodalGravity
    variable = ue_y
    boundary = pull
    gravity_value = -1
    mass = 1
  []
[]

[Materials]
  [elasticity_tensor_e]
    type = ComputeVariableIsotropicElasticityTensor
    youngs_modulus = E_phys
    poissons_ratio = poissons_ratio
    args = 'mat_den'
    base_name = "mech"
  []
  [elasticity_tensor_t]
    type = ComputeVariableIsotropicElasticityTensor
    youngs_modulus = E_phys
    poissons_ratio = poissons_ratio
    args = 'mat_den'
    base_name = "thermo"
  []
  
  [E_phys]
    type = DerivativeParsedMaterial
    # Emin + (density^penal) * (E0 - Emin)
    expression = '${Emin} + (mat_den ^ ${power}) * (${E0}-${Emin})'
    coupled_variables = 'mat_den'
    property_name = E_phys
  []
  [poissons_ratio]
    type = GenericConstantMaterial
    prop_names = poissons_ratio
    prop_values = 0.3
  []
  
  [stress_e] # calculate stress
    type = ComputeLinearElasticStress
    base_name = "mech"
  []
  [strain_e] # calculate strain
    type = ComputeSmallStrain
    base_name = "mech"
    displacements = 'ue_x ue_y ue_z'
  []
  
  [stress_t] # calculate stress
    type = ComputeLinearElasticStress
    base_name = "thermo"
  []
  [strain_t] # calculate strain
    type = ComputeSmallStrain
    base_name = "thermo"
    displacements = 'ut_x ut_y ut_z'
  []
  [./eigenstrain]
    type = ComputeEigenstrain
    #eigen_base = '0'
    eigen_base = '-1e-6 -1e-6 1e-6 0 0 0'
    eigenstrain_name = eigenstrain
  [../]
  
  [dc]
    type = ComplianceSensitivity
    design_density = mat_den
    youngs_modulus = E_phys
    incremental = false
    base_name = "mech"
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[UserObjects]
  [update]
    type = DensityUpdate
    density_sensitivity = Dc_elem
    design_density = mat_den
    volume_fraction = ${vol_frac}
    execute_on = TIMESTEP_BEGIN
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu'
  nl_abs_tol = 1e-10
  nl_max_its = 10
  line_search = none
  dt = 1.0
  num_steps = 10
[]

[Outputs]
  [out]
    type = Exodus
    time_step_interval = 10
  []
[]

