# H-field sub-app for staggered H-J problem
# Solves for magnetic field H using MFEMHypreAMS preconditioner (suitable for H(curl) problems)

pi = 6.2832
I = 96.0
E0 = 0.0001
Jc0 = 2.66e10
n = 21
mu0 = 1.25663706212e-6
r_outer = 0.009
H0 = 81805.6407047012
k = 0.6288

[Mesh]
  type = MFEMFileMesh
  file = 'multilayer100um.msh'
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
  [./H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
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
[]

# J is received from J_field_app via transfer - defined as AuxVariable
[AuxVariables]
  [./J]
    type = MFEMVariable
    fespace = HDivFESpace # try keep J on HDivFESpace
  [../]
  [./H_x]
    type = MFEMVariable
    fespace = H1FESpace
  [../]
  [./H_y]
    type = MFEMVariable
    fespace = H1FESpace
  [../]
  [./H_z]
    type = MFEMVariable
    fespace = H1FESpace
  [../]
  [./resistivity_out]
    type = MFEMVariable
    fespace = H1FESpace
  [../]
[]

[Functions]
  # Vector functions for tangential BC on H (HCurl space)
  [./outer_surface_H]
    type = ParsedVectorFunction
    expression_x = 't * (-y / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
    expression_y = 't * (x / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
    expression_z = '0.0'
  [../]
  # Resistivity functions for each material block (uses J_mag from transferred J)
  [./resistivity_tape]
    type = MFEMParsedFunction
    expression = '(${E0} / (${Jc0} * ${mu0})) * (J_mag / ${Jc0})^(${n}-1) * ((1.0 + sqrt(${k}*(H_x^2 + H_z^2) + H_y^2)/${H0})^${n})'
    symbol_names = 'J_mag H_x H_y H_z'
    symbol_values = 'J_mag H_x H_y H_z'
  [../]
  [./resistivity_copper]
    type = MFEMParsedFunction
    expression = '8e-9 / ${mu0}'
  [../]
  [./resistivity_silver]
    type = MFEMParsedFunction
    expression = '1e-8 / ${mu0}'
  [../]
  [./resistivity_substrate]
    type = MFEMParsedFunction
    expression = '1e-6 / ${mu0}'
  [../]
  [./resistivity_outer]
    type = MFEMParsedFunction
    expression = '1.0 / ${mu0}'
  [../]
  # Derivative term for nonlinear resistivity in tape: (n-1) * resistivity
  [./j_dresistivity_dj_tape]
    type = MFEMParsedFunction
    expression = '(${n}-1) * resistivity'
    symbol_names = 'resistivity'
    symbol_values = 'resistivity'
  [../]
  # Zero derivative for linear materials
  [./j_dresistivity_dj_zero]
    type = MFEMParsedFunction
    expression = '0.0'
  [../]
[]

[FunctorMaterials]
  [./resistivity_tape_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity'
    prop_values = 'resistivity_tape'
    block = 'tape'
  [../]
  [./resistivity_copper_bottom_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity'
    prop_values = 'resistivity_copper'
    block = 'copper_bottom'
  [../]
  [./resistivity_copper_top_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity'
    prop_values = 'resistivity_copper'
    block = 'copper_top'
  [../]
  [./resistivity_silver_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity'
    prop_values = 'resistivity_silver'
    block = 'silver'
  [../]
  [./resistivity_substrate_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity'
    prop_values = 'resistivity_substrate'
    block = 'substrate'
  [../]
  [./resistivity_outer_circle_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'resistivity'
    prop_values = 'resistivity_outer'
    block = 'outer_circle'
  [../]
  # Derivative materials for nonlinear resistivity (only needed for tape)
  [./j_dresistivity_dj_tape_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'j_dresistivity_dj'
    prop_values = 'j_dresistivity_dj_tape'
    block = 'tape'
  [../]
  [./j_dresistivity_dj_copper_bottom_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'j_dresistivity_dj'
    prop_values = 'j_dresistivity_dj_zero'
    block = 'copper_bottom'
  [../]
  [./j_dresistivity_dj_copper_top_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'j_dresistivity_dj'
    prop_values = 'j_dresistivity_dj_zero'
    block = 'copper_top'
  [../]
  [./j_dresistivity_dj_silver_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'j_dresistivity_dj'
    prop_values = 'j_dresistivity_dj_zero'
    block = 'silver'
  [../]
  [./j_dresistivity_dj_substrate_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'j_dresistivity_dj'
    prop_values = 'j_dresistivity_dj_zero'
    block = 'substrate'
  [../]
  [./j_dresistivity_dj_outer_circle_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'j_dresistivity_dj'
    prop_values = 'j_dresistivity_dj_zero'
    block = 'outer_circle'
  [../]
[]

[BCs]
  # H field tangential BCs (HCurl space) - MFEMVectorTangentialDirichletBC
  [./outer_surface_H]
    type = MFEMVectorTangentialDirichletBC
    variable = H
    vector_coefficient = outer_surface_H
    boundary = 'outer_circle_external'
  [../]
[]

[Kernels]
  # Faraday law: resistivity * curl H * curl test + dH/dt * test = 0
  # Using lagged nonlinearity with k_coefficient = resistivity(J_mag)
  # For tape block, curlu_dk_dcurlu_coefficient = j_dresistivity_dj
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
[]

[AuxKernels]
  # H components for output
  [./H_x]
    type = MFEMVectorComponentAux
    variable = H_x
    vector_coefficient = H
    component = 0
  [../]
  [./H_y]
    type = MFEMVectorComponentAux
    variable = H_y
    vector_coefficient = H
    component = 1
  [../]
  [./H_z]
    type = MFEMVectorComponentAux
    variable = H_z
    vector_coefficient = H
    component = 2
  [../]
  [./resistivity]
    type = MFEMScalarProjectionAux
    variable = resistivity_out
    coefficient = resistivity
  [../]
[]

# MultiApp definitions for the two sub-apps
[MultiApps]
  [./J_field_app]
    type = TransientMultiApp
    execute_on = 'TIMESTEP_END'
    input_files = 'tape_J_field.i'
    sub_cycling = true
  [../]
[]

# Transfers between sub-apps
[Transfers]
  # Transfer J from J_field_app to H_field_app (for resistivity calculation)
  [./J_to_H_app]
    type = MultiAppMFEMShapeEvaluationTransfer
    source_variables = 'J'
    variables = 'J'
    from_multi_app = J_field_app
    execute_on = 'TIMESTEP_END'
  [../]

  # Transfer H from H_field_app to J_field_app (for curl H in Ampere law)
  [./H_to_J_app]
    type = MultiAppMFEMShapeEvaluationTransfer
    source_variables = 'H'
    variables = 'H'
    to_multi_app = J_field_app
    execute_on = 'TIMESTEP_END'
  [../]
[]

# Preconditioner: MFEMHypreAMS for H(curl) problems
[Solvers]
  [./ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
  [../]
  [./main]
    type = MFEMHyprePCG
    preconditioner = ams
    l_tol = 1e-10
    l_max_its = 5000
    print_level = 1
  [../]
  [./newton]
    type = MFEMNewtonNonlinearSolver
    use_initial_guess = true
    max_its = 100
    abs_tol = 1e-9
    rel_tol = 1e-5
    print_level = 1
  [../]
[]

[Executioner]
  type = MFEMTransient
  start_time = 0.0
  end_time = 1.11
  dt = 0.0025
[]

[Outputs]
  [./out]
    type = MFEMParaViewDataCollection
    file_base = OutputData/tape_H_field
    vtk_format = ASCII
    time_step_interval = 40
  [../]
[]
