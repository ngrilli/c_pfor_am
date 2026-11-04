pi = 6.2832
I = 96.0
#r_inner = 1.0
r_outer = 2.0
E0 = 0.0001
Jc0 = 2.66e10
n = 21
mu0 = 1.25663706212e-6
H0 = 81805.6407047012
k = 0.6288

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

  [./lambda_tape_inlet]
    family = SCALAR
    order = FIRST
  [../]
  [./lambda_tape_outlet]
    family = SCALAR
    order = FIRST
  [../]
[]

[AuxVariables]
  [./resistivity]
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
  [./resistivity]
    type = MaterialRealAux
    variable = resistivity
    property = resistivity
  [../]
[]

[Materials]
  [./resistivity_inner]
    type = ParsedMaterial
    property_name = resistivity
    coupled_variables = 'Jx Jy Jz Hx Hy Hz'
    expression = '(${E0} / (${Jc0} * ${mu0})) * ((sqrt(Jx^2 + Jy^2 + Jz^2) / ${Jc0})^(${n}-1)) * ((1.0 + sqrt(${k}*(Hx^2 + Hz^2) + Hy^2)/${H0})^${n})'
    block = 'inner_block'
  [../]
  [./resistivity_outer]
    type = ParsedMaterial
    property_name = resistivity
    expression = '1.0 / ${mu0}'
    block = 'outer_block'
  [../]
[]

[BCs]
  [./sideHz]
    type = DirichletBC
    variable = Hz
    boundary = 'outer_circle'
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

#  [./back_front_Hx]
#    type = FunctionDirichletBC
#    variable = Hx
#    boundary = 'outer_back inner_back outer_front inner_front'
#    function = back_front_Hx
#  [../]
#  [./back_front_Hy]
#    type = FunctionDirichletBC
#    variable = Hy
#    boundary = 'outer_back inner_back outer_front inner_front'
#    function = back_front_Hy
#  [../]

  [./J_inlet]
    type = BoundaryIntegralValueConstraint
    variable = Jz
    lambda = lambda_tape_inlet
    phi0 = total_current
    boundary = 'inner_back'
  [../]
  [./J_outlet]
    type = BoundaryIntegralValueConstraint
    variable = Jz
    lambda = lambda_tape_outlet
    phi0 = total_current
    boundary = 'inner_front'
  [../]

#  [./J_zero]
#    type = FunctionDirichletBC
#    variable = Jz
#    function = '0.0'
#    boundary = 'outer_back outer_front'
#  [../]
[]

[Functions]
  [./outer_surface_Hx]
    type = ParsedFunction
    expression = 'min(t,1.0) * (-y / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
  [../]
  [./outer_surface_Hy]
    type = ParsedFunction
    expression = 'min(t,1.0) * (x / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
  [../]
#  [./back_front_Hx]
#    type = ParsedFunction
#    expression = 'if(sqrt(x^2 + y^2)<${r_inner},(-y / sqrt(x^2 + y^2)) * ${I} * sqrt(x^2 + y^2) / (2.0 * ${pi} * ${r_inner}^2),(-y / sqrt(x^2 + y^2)) * ${I} / (2.0 * ${pi} * sqrt(x^2 + y^2)))'
#  [../]
#  [./back_front_Hy]
#    type = ParsedFunction
#    expression = 'if(sqrt(x^2 + y^2)<${r_inner},(x / sqrt(x^2 + y^2)) * ${I} * sqrt(x^2 + y^2) / (2.0 * ${pi} * ${r_inner}^2),(x / sqrt(x^2 + y^2)) * ${I} / (2.0 * ${pi} * sqrt(x^2 + y^2)))'
#  [../]
[]

[Postprocessors]
  [./time]
    type = TimePostprocessor
    execute_on = 'timestep_begin'
  [../]
  [./total_current]
    type = ParsedPostprocessor
    expression = '16000000000.0 * min(time,1.0)'
    pp_names = 'time'
    execute_on = 'timestep_begin'
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
  
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  line_search = 'none'
  
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8

  nl_max_its = 30
  l_max_its = 200
  
  start_time = 0.0
  end_time = 2.0
  dt = 0.005
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]
