pi = 6.2832
I = 1.0
r_inner = 1.0
r_outer = 2.0

[Mesh]
  [./cylinder]
    type = FileMeshGenerator
    file = 'concentric_cylinder_2D.msh'
  [../]
[]

[Variables]
  [H]
    family = NEDELEC_ONE
    order = FIRST
  []
[]

[Kernels]
  [curlcurl]
    type = CurlCurlField
    variable = H
  []
  [./TimeDerivative]
    type = VectorTimeDerivative
    variable = H
  [../]
[]

[BCs]
    # Outer cylinder boundary - tangential H field
    [./outer_surface]
      type = VectorCurlPenaltyDirichletBC
      variable = H
      function_x = outer_surface_Hx
      function_y = outer_surface_Hy
      penalty = 1e8
      boundary = outer_circle
    [../]

    # Inner cylinder boundary - tangential H field
    [./inner_surface]
      type = VectorCurlPenaltyDirichletBC
      variable = H
      function_x = inner_surface_Hx
      function_y = inner_surface_Hy
      penalty = 1e8
      boundary = inner_circle
    [../]

    # Back/front boundaries - tangential H field
    #[./back_front]
    #  type = VectorCurlPenaltyDirichletBC
    #  variable = H
    #  function_x = back_front_Hx
    #  function_y = back_front_Hy
    #  penalty = 1e8
    #  boundary = 'back front'
    #[../]
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
  [SMP]
    type = SMP
    full = true
  []
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