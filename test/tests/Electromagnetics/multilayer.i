pi = 6.2832
I = 96.0
tape_surface = 18.0e-9
E0 = 0.0001
Jc0 = 2.58e10
n = 21
mu0 = 1.25663706212e-6
r_outer = 0.009

[Mesh]
  [./cylinder]
    type = FileMeshGenerator
    file = 'multilayer.msh'
  [../]
[]

[Variables]
  [./Hx]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Hy]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Hz]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Jx]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Jy]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Jz]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
[]

[Kernels]
  [./curlcurl_x]
    type = CurlCurlComponent
    variable = Hx
    H1 = Hz
    H2 = Hy
    component = 0
    resistivity = resistivity
  [../]
  [./curlcurl_y]
    type = CurlCurlComponent
    variable = Hy
    H1 = Hx
    H2 = Hz
    component = 1
    resistivity = resistivity
  [../]
  [./curlcurl_z]
    type = CurlCurlComponent
    variable = Hz
    H1 = Hy
    H2 = Hx
    component = 2
    resistivity = resistivity
  [../]

  [./TimeDerivative_x]
    type = TimeDerivative
    variable = Hx
  [../]
  [./TimeDerivative_y]
    type = TimeDerivative
    variable = Hy
  [../]
  [./TimeDerivative_z]
    type = TimeDerivative
    variable = Hz
  [../]

  [./curl_x]
    type = CoupledCurlComponent
    variable = Jx
    Hx = Hx
    Hy = Hy
    Hz = Hz
    component = 0
  [../]
  [./curl_y]
    type = CoupledCurlComponent
    variable = Jy
    Hx = Hx
    Hy = Hy
    Hz = Hz
    component = 1
  [../]
  [./curl_z]
    type = CoupledCurlComponent
    variable = Jz
    Hx = Hx
    Hy = Hy
    Hz = Hz
    component = 2
  [../]

  [./Jx]
    type = Reaction
    variable = Jx
    rate = -1.0
  [../]
  [./Jy]
    type = Reaction
    variable = Jy
    rate = -1.0
  [../]
  [./Jz]
    type = Reaction
    variable = Jz
    rate = -1.0
  [../]
[]

[AuxKernels]
[]

[Materials]
  # resistivity is divided by magnetic permeability
  [./resistivity_tape]
    type = ParsedMaterial
    property_name = resistivity
    coupled_variables = 'Jx Jy Jz'
    expression = '(${E0} / (${Jc0} * ${mu0})) * (sqrt(Jx^2 + Jy^2 + Jz^2) / ${Jc0})^(${n}-1)'
    block = 'silver' #'tape'
  [../]
  [./resistivity_copper]
    type = ParsedMaterial
    property_name = resistivity
    expression = '8e-9 / ${mu0}'
    block = 'copper_bottom copper_top'
  [../]
#  [./resistivity_silver]
#    type = ParsedMaterial
#    property_name = resistivity
#    expression = '1e-8 / ${mu0}'
#    block = 'silver'
#  [../]
  [./resistivity_substrate]
    type = ParsedMaterial
    property_name = resistivity
    expression = '1e-6 / ${mu0}'
    block = 'substrate'
  [../]
  [./resistivity_outer]
    type = ParsedMaterial
    property_name = resistivity
    expression = '1.0 / ${mu0}'
    block = 'outer_circle'
  [../]
[]

[BCs]
  [./outer_surface_Hx]
    type = FunctionDirichletBC
    variable = Hx
    boundary = outer_circle_external
    function = outer_surface_Hx
  [../]
  [./outer_surface_Hy]
    type = FunctionDirichletBC
    variable = Hy
    boundary = outer_circle_external
    function = outer_surface_Hy
  [../]

  [./J_inlet_outlet]
    type = FunctionDirichletBC
    variable = Jz
    function = 't * ${I} / ${tape_surface}'
    boundary = 'silver_inlet silver_outlet' #'tape_inlet tape_outlet'
  [../]
  [./Jz_zero]
    type = FunctionDirichletBC
    variable = Jz
    function = '0.0'
    boundary = 'outer_circle_inlet outer_circle_outlet copper_top_inlet copper_top_outlet copper_bottom_inlet copper_bottom_outlet substrate_inlet substrate_outlet'
  [../]
  [./Jx_zero]
    type = FunctionDirichletBC
    variable = Jx
    function = '0.0'
    boundary = 'silver_inlet silver_outlet outer_circle_inlet outer_circle_outlet copper_top_inlet copper_top_outlet copper_bottom_inlet copper_bottom_outlet substrate_inlet substrate_outlet'
  [../]
  [./Jy_zero]
    type = FunctionDirichletBC
    variable = Jy
    function = '0.0'
    boundary = 'silver_inlet silver_outlet outer_circle_inlet outer_circle_outlet copper_top_inlet copper_top_outlet copper_bottom_inlet copper_bottom_outlet substrate_inlet substrate_outlet'
  [../]
[]

[Functions]
  [./outer_surface_Hx]
    type = ParsedFunction
    expression = 't * (-y / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
  [../]
  [./outer_surface_Hy]
    type = ParsedFunction
    expression = 't * (x / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg          51'
  line_search = 'none'
  
  nl_rel_tol = 1e-3
  nl_abs_tol = 1e-3

  nl_max_its = 50
  l_max_its = 100
  
  start_time = 0.0
  end_time = 1.0
  dt = 0.001
[]

[Outputs]
  [./out]
    type = Exodus
    time_step_interval = 10
  [../]
[]
