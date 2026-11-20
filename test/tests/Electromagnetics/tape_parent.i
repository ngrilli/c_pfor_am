# Parent app for staggered H-J problem using MultiApp
# Coordinates H-field and J-field sub-apps with transfers between them

pi = 6.2832
I = 96.0
tape_surface = 6.0e-9
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
  solve = false
[]

# Define FESpaces needed for transfers (same as sub-apps)
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
  [./H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  [../]
[]

# AuxVariables for receiving transferred data from sub-apps
[AuxVariables]
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
  [./J_x]
    type = MFEMVariable
    fespace = H1FESpace
  [../]
  [./J_y]
    type = MFEMVariable
    fespace = H1FESpace
  [../]
  [./J_z]
    type = MFEMVariable
    fespace = H1FESpace
  [../]
  [./resistivity_out]
    type = MFEMVariable
    fespace = H1FESpace
  [../]
[]

# Functions for boundary conditions (shared with sub-apps)
[Functions]
  # Vector functions for tangential BC on H (HCurl space)
  [./outer_surface_H]
    type = ParsedVectorFunction
    expression_x = 't * (-y / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
    expression_y = 't * (x / ${r_outer}) * ${I} / (2.0 * ${pi} * ${r_outer})'
    expression_z = '0.0'
  [../]
  # Vector functions for normal BC on J (HDiv space)
  [./J_zero]
    type = ParsedVectorFunction
    expression_x = '0.0'
    expression_y = '0.0'
    expression_z = '0.0'
  [../]
  [./J_inlet]
    type = ParsedVectorFunction
    expression_x = '0.0'
    expression_y = '0.0'
    expression_z = 't * ${I} / ${tape_surface}'
  [../]
  # Resistivity functions for each material block
  [./resistivity_tape]
    type = MFEMParsedFunction
    expression = '(${E0} / (${Jc0} * ${mu0})) * (J_mag / ${Jc0})^(${n}-1) * ((1.0 + sqrt(${k}*(H_x^2 + H_z^2) + H_y^2)/${H0})^${n})'
    symbol_names = 'J_mag H_mag H_x H_y H_z'
    symbol_values = 'J_mag H_mag H_x H_y H_z'
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

# MultiApp definitions for the two sub-apps
[MultiApps]
  [./H_field_app]
    type = FullSolveMultiApp
    input_files = 'tape_H_field.i'
    execute_on = 'TIMESTEP_BEGIN'
    sub_cycling = true
  [../]
  [./J_field_app]
    type = FullSolveMultiApp
    input_files = 'tape_J_field.i'
    execute_on = 'TIMESTEP_BEGIN'
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
    to_multi_app = H_field_app
    execute_on = 'TIMESTEP_END'
  [../]

  # Transfer H from H_field_app to J_field_app (for curl H in Ampere law)
  [./H_to_J_app]
    type = MultiAppMFEMShapeEvaluationTransfer
    source_variables = 'H'
    variables = 'H'
    from_multi_app = H_field_app
    to_multi_app = J_field_app
    execute_on = 'TIMESTEP_END'
  [../]

  # Also transfer to parent for output/postprocessing
  [./H_from_H_app]
    type = MultiAppMFEMShapeEvaluationTransfer
    source_variables = 'H'
    variables = 'H'
    from_multi_app = H_field_app
    execute_on = 'TIMESTEP_END'
  [../]

  [./J_from_J_app]
    type = MultiAppMFEMShapeEvaluationTransfer
    source_variables = 'J'
    variables = 'J'
    from_multi_app = J_field_app
    execute_on = 'TIMESTEP_END'
  [../]
[]

# Transient executioner (same as sub-apps)
[Executioner]
  type = MFEMTransient
  start_time = 0.0
  end_time = 1.0
  dt = 0.002
[]

[Outputs]
  [./out]
    type = MFEMParaViewDataCollection
    file_base = OutputData/tape_parent
    vtk_format = ASCII
    time_step_interval = 10
  [../]
[]