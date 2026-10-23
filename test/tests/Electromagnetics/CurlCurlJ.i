pi = 6.2832
I = 1.0
r_inner = 1.0
r_outer = 2.0
E0 = 0.0001
Jc0 = 2.58e10
n = 21

[Mesh]
  [./cylinder]
    type = FileMeshGenerator
    file = 'cylinder_3D.msh'
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
  [./aux_Jx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./aux_Jy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./aux_Jz]
    order = CONSTANT
    family = MONOMIAL
  [../]
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
  [./aux_Jx]
    type = CurrentDensityComponent
    variable = aux_Jx
    Hx = Hx
    Hy = Hy
    Hz = Hz
    component = 0
  [../]
  [./aux_Jy]
    type = CurrentDensityComponent
    variable = aux_Jy
    Hx = Hx
    Hy = Hy
    Hz = Hz
    component = 1
  [../]
  [./aux_Jz]
    type = CurrentDensityComponent
    variable = aux_Jz
    Hx = Hx
    Hy = Hy
    Hz = Hz
    component = 2
  [../]
[]

[Materials]
  [./resistivity_inner]
    type = ParsedMaterial
    property_name = resistivity
    coupled_variables = 'Jx Jy Jz'
    expression = '(${E0} / ${Jc0}) * (sqrt(Jx^2 + Jy^2 + Jz^2) / ${Jc0})^(${n}-1)'
    block = 'inner_block'
  [../]
  [./resistivity_outer]
    type = ParsedMaterial
    property_name = resistivity
    coupled_variables = 'Jx Jy Jz'
    expression = '1.0'
    block = 'outer_block'
  [../]
[]

[BCs]
  [./sideHz]
    type = DirichletBC
    variable = Hz
    boundary = 'outer_circle outer_back inner_back outer_front inner_front'
    value = 0.0
  [../]

  [./outer_surface_Hx]
    type = FunctionDirichletBC
    variable = Hx
    boundary = outer_circle
    function = outer_surface_Hx
  [../]
  [./outer_surface_Hy]
    type = FunctionDirichletBC
    variable = Hy
    boundary = outer_circle
    function = outer_surface_Hy
  [../]

  [./back_front_Hx]
    type = FunctionDirichletBC
    variable = Hx
    boundary = 'outer_back inner_back outer_front inner_front'
    function = back_front_Hx
  [../]
  [./back_front_Hy]
    type = FunctionDirichletBC
    variable = Hy
    boundary = 'outer_back inner_back outer_front inner_front'
    function = back_front_Hy
  [../]

  [./J_inlet]
    type = FunctionDirichletBC
    variable = Jz
    function = '${I} / (${pi} * ${r_inner}^2)'
    boundary = 'inner_back inner_front'
  [../]
  [./J_zero]
    type = FunctionDirichletBC
    variable = Jz
    function = '0.0'
    boundary = 'outer_back outer_front'
  [../]
[]

[Functions]
  [./outer_surface_Hx]
    type = ParsedFunction
    expression = '(-y / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
  [../]
  [./outer_surface_Hy]
    type = ParsedFunction
    expression = '(x / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
  [../]
  [./back_front_Hx]
    type = ParsedFunction
    expression = 'if(sqrt(x^2 + y^2)<${r_inner},(-y / sqrt(x^2 + y^2)) * ${I} * sqrt(x^2 + y^2) / (2.0 * ${pi} * ${r_inner}^2),(-y / sqrt(x^2 + y^2)) * ${I} / (2.0 * ${pi} * sqrt(x^2 + y^2)))'
  [../]
  [./back_front_Hy]
    type = ParsedFunction
    expression = 'if(sqrt(x^2 + y^2)<${r_inner},(x / sqrt(x^2 + y^2)) * ${I} * sqrt(x^2 + y^2) / (2.0 * ${pi} * ${r_inner}^2),(x / sqrt(x^2 + y^2)) * ${I} / (2.0 * ${pi} * sqrt(x^2 + y^2)))'
  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP # FDP
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
  
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6

  nl_max_its = 20
  l_max_its = 100
  
  start_time = 0.0
  end_time = 1.0
  dt = 0.005
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]
