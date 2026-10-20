pi = 6.2832
I = 1.0
r_inner = 1.0
r_outer = 2.0

[Mesh]
  [./cylinder]
    type = FileMeshGenerator
    file = 'concentric_cylinder_3D.msh'
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

  // TO DO: need a kernel similar to CurlComponent but with all three components of H as coupled variables
  [./CurlEComponent1_air]
    type = PowerLawCurlComponent
    variable = Jx
    J1 = Jz
    J2 = Jy
    E0 = 1.0
    J0 = 1.0
    n = 1.0
    component = 0
  [../]
  [./CurlEComponent2_air]
    type = PowerLawCurlComponent
    variable = Jy
    J1 = Jx
    J2 = Jz
    E0 = 1.0
    J0 = 1.0
    n = 1.0
    component = 1
  [../]
  [./CurlEComponent3_air]
    type = PowerLawCurlComponent
    variable = Jz
    J1 = Jy
    J2 = Jx
    E0 = 1.0
    J0 = 1.0
    n = 1.0
    component = 2
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
  [./resistivity]
    type = GenericConstantMaterial
    prop_names = 'resistivity'
    prop_values = '1.0'
  [../]
[]

[BCs]
  [./sideHz]
    type = DirichletBC
    variable = Hz
    boundary = 'outer_circle inner_circle back front'
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

  [./inner_surface_Hx]
    type = FunctionDirichletBC
    variable = Hx
    boundary = inner_circle
    function = inner_surface_Hx
  [../]
  [./inner_surface_Hy]
    type = FunctionDirichletBC
    variable = Hy
    boundary = inner_circle
    function = inner_surface_Hy
  [../]

  [./back_front_Hx]
    type = FunctionDirichletBC
    variable = Hx
    boundary = 'back front'
    function = back_front_Hx
  [../]
  [./back_front_Hy]
    type = FunctionDirichletBC
    variable = Hy
    boundary = 'back front'
    function = back_front_Hy
  [../]
[]

[Functions]
  [./outer_surface_Hx]
    type = ParsedFunction
    expression = '(-y / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
  [../]
  [./inner_surface_Hx]
    type = ParsedFunction
    expression = '(-y / ${r_inner}) * ${I} / (2.0 * ${pi} * ${r_inner})'
  [../]
  [./outer_surface_Hy]
    type = ParsedFunction
    expression = '(x / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
  [../]
  [./inner_surface_Hy]
    type = ParsedFunction
    expression = '(x / ${r_inner}) * ${I} / (2.0 * ${pi} * ${r_inner})'
  [../]
  [./back_front_Hx]
    type = ParsedFunction
    expression = '(-y / sqrt(x^2 + y^2)) * ${I} / (2.0 * ${pi} * sqrt(x^2 + y^2))'
  [../]
  [./back_front_Hy]
    type = ParsedFunction
    expression = '(x / sqrt(x^2 + y^2)) * ${I} / (2.0 * ${pi} * sqrt(x^2 + y^2))'
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
