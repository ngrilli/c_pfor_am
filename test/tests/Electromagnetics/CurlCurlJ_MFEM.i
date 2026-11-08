# MFEM version of CurlCurlJ.i
# Faraday law: resistivity * curl H * curl test + dH/dt * test = 0
# Ampere law: curl H test - J test = 0

pi = 6.2832
I = 96.0
r_inner = 1.0
r_outer = 2.0
E0 = 0.0001
Jc0 = 2.66e10
n = 21

[Mesh]
  type = MFEMMesh
  file = 'cylinder_MFEM.msh'
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [./HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  [../]
  [./HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = FIRST
  [../]
[]

[Variables]
  [./H]
    type = MFEMVariable
    fespace = HCurlFESpace
  [../]
  [./J]
    type = MFEMVariable
    fespace = HDivFESpace
  [../]
[]

[Functions]
  # Vector functions for tangential BC on H (HCurl space)
  [./outer_surface_H]
    type = ParsedVectorFunction
    expression_x = '(-y / ${r_outer}) * t * ${I} / (2.0 * ${pi} * ${r_outer})'
    expression_y = '(x / ${r_outer}) * t * ${I} / (2.0 * ${pi} * ${r_outer})'
    expression_z = '0.0'
  [../]
  [./back_front_H]
    type = ParsedVectorFunction
    expression_x = 'if(sqrt(x^2 + y^2)<${r_inner},(-y / sqrt(x^2 + y^2)) * ${I} * sqrt(x^2 + y^2) / (2.0 * ${pi} * ${r_inner}^2),(-y / sqrt(x^2 + y^2)) * ${I} / (2.0 * ${pi} * sqrt(x^2 + y^2)))'
    expression_y = 'if(sqrt(x^2 + y^2)<${r_inner},(x / sqrt(x^2 + y^2)) * ${I} * sqrt(x^2 + y^2) / (2.0 * ${pi} * ${r_inner}^2),(x / sqrt(x^2 + y^2)) * ${I} / (2.0 * ${pi} * sqrt(x^2 + y^2)))'
    expression_z = '0.0'
  [../]
  # Vector functions for normal BC on J (HDiv space)
  [./J_inner]
    type = ParsedVectorFunction
    expression_x = '0.0'
    expression_y = '0.0'
    expression_z = 't * ${I} / (${pi} * ${r_inner}^2)'
  [../]
  [./J_zero]
    type = ParsedVectorFunction
    expression_x = '0.0'
    expression_y = '0.0'
    expression_z = '0.0'
  [../]
  # Resistivity as function of |J| (J_mag is automatically available for HDiv variables)
  [./resistivity_inner]
    type = MFEMParsedFunction
    expression = '(${E0} / ${Jc0}) * (J_mag / ${Jc0})^(${n}-1)'
    symbol_names = 'J_mag'
    symbol_values = 'J_mag'
  [../]
  [./resistivity_outer]
    type = MFEMParsedFunction
    expression = '1.0'
  [../]
  # Derivative term for nonlinear resistivity: (n-1) * resistivity
  [./j_dresistivity_dj]
    type = MFEMParsedFunction
    expression = '(${n}-1) * resistivity'
    symbol_names = 'resistivity'
    symbol_values = 'resistivity'
  [../]
[]

[FunctorMaterials]
  [./resistivity_inner_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity'
    prop_values = 'resistivity_inner'
    block = 'inner_block'
  [../]
  [./resistivity_outer_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity'
    prop_values = 'resistivity_outer'
    block = 'outer_block'
  [../]
[]

[BCs]
  # H field tangential BCs (HCurl space) - MFEMVectorTangentialDirichletBC
  [./outer_surface_H]
    type = MFEMVectorTangentialDirichletBC
    variable = H
    vector_coefficient = outer_surface_H
    boundary = 'outer_circle'
  [../]
  #[./back_front_H]
  #  type = MFEMVectorTangentialDirichletBC
  #  variable = H
  #  vector_coefficient = back_front_H
  #  boundary = 'outer_back inner_back outer_front inner_front'
  #[../]

  # J field normal BCs (HDiv space) - MFEMVectorNormalDirichletBC
  [./J_inlet]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_inner
    boundary = 'inner_back inner_front'
  [../]
  [./J_zero]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_zero
    boundary = 'outer_back outer_front outer_circle'
  [../]
[]

[Kernels]
  # Faraday law: resistivity * curl H * curl test + dH/dt * test = 0
  # Note: MFEMNLCurlCurlKernel expects nonlinearity in |curl H|, but our resistivity depends on |J|.
  # Using lagged nonlinearity (curlu_dk_dcurlu_coefficient = 0) with k_coefficient = resistivity(J_mag).
  # For full Newton coupling, a custom kernel would be needed.
  [./curlcurl]
    type = MFEMNLCurlCurlKernel
    variable = H
    k_coefficient = resistivity
    curlu_dk_dcurlu_coefficient = j_dresistivity_dj
  [../]
  [./dHdt]
    type = MFEMTimeDerivativeVectorFEMassKernel
    variable = H
    coefficient = 1.0
  [../]

  # Ampere law: curl H test - J test = 0
  # MFEMVectorFECurlKernel: trial=H (HCurl), test=HDiv (curl H in HDiv tested with HDiv test function)
  # MFEMVectorFEMassKernel: trial=J (HDiv), test=HDiv, coefficient=-1
  [./curlH]
    type = MFEMVectorFECurlKernel
    variable = H
    coefficient = 1.0
  [../]
  [./J_term]
    type = MFEMVectorFEMassKernel
    variable = J
    coefficient = -1.0
  [../]
[]

[Preconditioner]
  [./ams_H]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
  [../]
  [./ams_J]
    type = MFEMHypreAMS
    fespace = HDivFESpace
  [../]
[]

[Solvers]
  [./mumps]
    type = MFEMMUMPS
    print_level = 0
  [../]
  #[./pcg_H]
  #  type = MFEMHyprePCG
  #  preconditioner = ams_H
  #  l_tol = 1e-12
  #  l_max_its = 100
  #  print_level = 0
  #[../]
  #[./pcg_J]
  #  type = MFEMHyprePCG
  #  preconditioner = ams_J
  #  l_tol = 1e-12
  #  l_max_its = 100
  #  print_level = 0
  #[../]
  [./newton]
    type = MFEMNewtonNonlinearSolver
    max_its = 200
    abs_tol = 1e-6
    rel_tol = 1e-6
    print_level = 1
  [../]
[]

[Executioner]
  type = MFEMTransient
  start_time = 0.0
  end_time = 1.0
  dt = 0.005
[]

[Outputs]
  [./out]
    type = MFEMParaViewDataCollection
    file_base = OutputData/CurlCurlJ_MFEM
    vtk_format = ASCII
  [../]
[]