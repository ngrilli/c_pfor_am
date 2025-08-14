# example with two sensitivities

vol_frac = 0.4
E0 = 1e5
Emin = 1e-4
power = 2

[GlobalParams]
  #displacements = 'disp_x disp_y disp_z'
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
    coord = '0 0 5'
  []
[]

[Variables]
  [Dc]
    initial_condition = -1.0
  []
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
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
  [./disp_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
    base_name = "mech"
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./disp_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    base_name = "mech"
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./disp_z]
    type = StressDivergenceTensors
    variable = disp_z
    component = 2
    base_name = "mech"
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[BCs]
  [no_x]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0.0
  []
  [no_y]
    type = DirichletBC
    variable = disp_y
    boundary = right
    value = 0.0
  []
  [no_z]
    type = DirichletBC
    variable = disp_z
    boundary = right
    value = 0.0
  []
  [boundary_penalty]
    type = ADRobinBC
    variable = Dc
    boundary = 'left top front back'
    coefficient = 10
  []
[]

[NodalKernels]
  [pull]
    type = NodalGravity
    variable = disp_y
    boundary = pull
    gravity_value = -1
    mass = 1
    #base_name = "mech"
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeVariableIsotropicElasticityTensor
    youngs_modulus = E_phys
    poissons_ratio = poissons_ratio
    args = 'mat_den'
    base_name = "mech"
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
  
  [stress] # calculate stress
    type = ComputeLinearElasticStress
    base_name = "mech"
  []
  [strain] # calculate strain
    type = ComputeSmallStrain
    base_name = "mech"
    displacements = 'disp_x disp_y disp_z'
  []
  
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

