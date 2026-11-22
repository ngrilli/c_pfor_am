# J-field sub-app for staggered H-J problem
# Solves for current density J using MFEMHypreADS preconditioner (suitable for H(div) problems)

I = 96.0
tape_surface = 6.0e-9

[Mesh]
  type = MFEMFileMesh
  file = 'multilayer100um.msh'
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
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
  [./HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  [../]
[]

[Variables]
  [./J]
    type = MFEMVariable
    fespace = HDivFESpace
  [../]
[]

# H is received from H_field_app via transfer - defined as AuxVariable
[AuxVariables]
  [./H]
    type = MFEMVariable
    fespace = HCurlFESpace  # H is kept as HCurl for MFEMVectorFECurlKernel
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
  # Vector functions for normal BC on J (HDiv space)
  [./J_inlet]
    type = ParsedVectorFunction
    expression_x = '0.0'
    expression_y = '0.0'
    expression_z = 't * ${I} / ${tape_surface}'
  [../]
[]

[BCs]
  # J field normal BCs (HDiv space) - MFEMVectorNormalDirichletBC
  # Jz is constrained via boundary integral on tape_inlet and tape_outlet
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
    type = MFEMVectorFEMassKernel # may have a -1 by default on left hand side of equation, see ex22p.cpp in MFEM
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

# Preconditioner: MFEMHypreADS for H(div) problems
[Solvers]
  inactive = 'jacobi main newton gmres'
  [./ads]
    type = MFEMHypreADS
    fespace = HDivFESpace
  [../]
  [./main]
    type = MFEMCGSolver
    preconditioner = jacobi
    l_tol = 1e-8
    l_max_its = 1000
    print_level = 2
  [../]
  [./HyprePCG]
    type = MFEMHyprePCG
    preconditioner = ads
    l_tol = 1e-10
    l_max_its = 10000
    print_level = 1
  [../]

  [./newton]
    type = MFEMNewtonNonlinearSolver
    use_initial_guess = true
    max_its = 1000
    abs_tol = 1e-12
    rel_tol = 1e-8
    print_level = 1
  [../]
  [./jacobi]
    type = MFEMOperatorJacobiSmoother
  [../]
  [./gmres]
    type = MFEMHypreGMRES
    preconditioner = boomeramg
    l_tol = 1e-8
    l_max_its = 1000
  [../]
[]

[Executioner]
  type = MFEMTransient
  start_time = 0.0
  end_time = 1.11
  dt = 0.0025
[]

[VectorPostprocessors]
  [./line_sample]
    type = MFEMVariableLineValueSampler
    variables = 'J_z'
    start_point = '-0.004 0.0000245 0.000300'
    end_point = '0.004 0.0000245 0.000300'
    num_points = 121
  [../]
[]

[Outputs]
  [./out]
    type = MFEMParaViewDataCollection
    file_base = OutputData/tape_J_field
    vtk_format = ASCII
    time_step_interval = 40
  [../]
  [./CSV]
    type = CSV
    execute_on = 'timestep_end'
    file_base = OutputData/tape_line_plot
    time_step_interval = 40
  [../]
[]
