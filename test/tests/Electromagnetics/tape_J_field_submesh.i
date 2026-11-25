# J-field sub-app for Tape geometry
# Solves for current density J on submesh (tape + copper_bottom + copper_top + silver + substrate blocks) using MFEMHypreADS preconditioner (suitable for H(div) problems)
# H field is received from H_field_app via MultiApp transfer
# Zero normal BC on copper_bottom_sides and copper_top_sides boundaries

I = 96.0
tape_surface = 6.0e-9

[Mesh]
  type = MFEMFileMesh
  file = 'multilayer100um.msh'
[]

[Problem]
  type = MFEMProblem
[]

# Submesh for tape, copper_bottom, copper_top, silver, and substrate blocks where J is defined
[SubMeshes]
  [./conductor_submesh]
    type = MFEMDomainSubMesh
    block = 'tape copper_bottom copper_top silver substrate'
    submesh_boundary = conductor_surface
  [../]
[]

[FESpaces]
  # HDivFESpace (RT) on submesh for J variable
  [./HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = FIRST
    submesh = conductor_submesh
  [../]
  # HCurlFESpace (ND) on submesh for H aux variable (received from main app)
  [./HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    submesh = conductor_submesh
  [../]
  # H1FESpace on submesh for J component output
  [./H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = conductor_submesh
  [../]
[]

[Variables]
  [./J]
    type = MFEMVariable
    fespace = HDivFESpace
  [../]
[]

# H is received from H_field_app via transfer - defined as AuxVariable on submesh
[AuxVariables]
  [./H]
    type = MFEMVariable
    fespace = HCurlFESpace  # H kept as HCurl for MFEMVectorFECurlKernel
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
[]

[Functions]
  # Zero vector function for normal BC on J (HDiv space) on copper_bottom_sides and copper_top_sides
  [./J_zero_normal]
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
[]

[BCs]
  # J field normal BCs (HDiv space) - MFEMVectorNormalDirichletBC
  # Zero normal BC on copper_bottom sides and copper_top sides boundaries
  [./copper_bottom_sides_J_zero]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_zero_normal
    boundary = 'copper_bottom_right copper_bottom_left copper_bottom_bottom'
  [../]
  [./copper_top_sides_J_zero]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_zero_normal
    boundary = 'copper_top_right copper_top_left copper_top_top'
  [../]
  [./silver_sides_J_zero]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_zero_normal
    boundary = 'silver_right silver_left'
  [../]
  [./substrate_sides_J_zero]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_zero_normal
    boundary = 'substrate_right substrate_left'
  [../]
  [./tape_sides_J_zero]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_zero_normal
    boundary = 'tape_right tape_left'
  [../]
  [./J_tape_inlet]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_inlet
    boundary = 'tape_inlet'
  [../]
  [./J_tape_outlet]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_inlet
    boundary = 'tape_outlet'
  [../]
[]

[Kernels]
  # Ampere law: curl H test - J test = 0
  # MFEMVectorFECurlKernel: trial=H (HCurl), test=HDiv (curl H in HDiv tested with HDiv test function)
  # MFEMVectorFEMassKernel: trial=J (HDiv), test=HDiv
  [./curlH]
    type = MFEMVectorFECurlKernel
    trial_variable = H
    variable = J
    coefficient = -1.0
  [../]
  [./J_term]
    type = MFEMVectorFEMassKernel
    variable = J
    coefficient = 1.0
  [../]
[]

[AuxKernels]
  # J components for output
  [./J_x]
    type = MFEMVectorComponentAux
    variable = J_x
    vector_coefficient = J
    component = 0
  [../]
  [./J_y]
    type = MFEMVectorComponentAux
    variable = J_y
    vector_coefficient = J
    component = 1
  [../]
  [./J_z]
    type = MFEMVectorComponentAux
    variable = J_z
    vector_coefficient = J
    component = 2
  [../]
[]

# Preconditioner: MFEMHypreADS for H(div) problems (matching cable_J_field.i settings)
[Solvers]
  [./ads]
    type = MFEMHypreADS
    fespace = HDivFESpace
  [../]
  [./HyprePCG]
    type = MFEMHyprePCG
    preconditioner = ads
    l_tol = 1e-10
    l_max_its = 10000
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
    file_base = OutputData/tape_J_field_submesh
    vtk_format = ASCII
    time_step_interval = 40
    submesh = conductor_submesh
  [../]
[]