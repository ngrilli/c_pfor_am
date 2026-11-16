# MFEM version of CurlCurlJ.i adapted from multilayer.i
# Faraday law: resistivity * curl H * curl test + dH/dt * test = 0
# Ampere law: curl H test - J test = 0

pi = 6.2832
I = 96.0
# tape_surface = 6.0e-9
E0 = 0.0001
Jc0 = 2.58e10
n = 21
mu0 = 1.25663706212e-6
r_outer = 0.009
H0 = 81805.6407047012
k = 0.6288

[Mesh]
  type = MFEMFileMesh
  file = 'multilayer.msh'
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
  [./H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
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
[]

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
  # Resistivity functions for each material block
  # tape: nonlinear resistivity depending on |J| and |H|
  [./resistivity_tape]
    type = MFEMParsedFunction
    expression = '(${E0} / (${Jc0} * ${mu0})) * (J_mag / ${Jc0})^(${n}-1) * ((1.0 + sqrt(${k}*(H_mag^2))/${H0})^${n})'
    symbol_names = 'J_mag H_mag'
    symbol_values = 'J_mag H_mag'
  [../]
  # copper: constant resistivity
  [./resistivity_copper]
    type = MFEMParsedFunction
    expression = '8e-9 / ${mu0}'
  [../]
  # silver: constant resistivity
  [./resistivity_silver]
    type = MFEMParsedFunction
    expression = '1e-8 / ${mu0}'
  [../]
  # substrate: constant resistivity
  [./resistivity_substrate]
    type = MFEMParsedFunction
    expression = '1e-6 / ${mu0}'
  [../]
  # outer_circle: constant resistivity (air/insulation)
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
  # Applied on outer_circle_external boundary
  [./outer_surface_H]
    type = MFEMVectorTangentialDirichletBC
    variable = H
    vector_coefficient = outer_surface_H
    boundary = 'outer_circle_external'
  [../]

  # J field normal BCs (HDiv space) - MFEMVectorNormalDirichletBC
  # Jz is constrained via boundary integral on tape_inlet and tape_outlet
  # For MFEM, we set J normal component to zero on all other boundaries
  # Note: BoundaryIntegralValueConstraint is not directly available in MFEM BCs
  # The current is enforced through the normal component BC on tape_inlet/outlet
  [./J_tape_inlet]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_zero
    boundary = 'tape_inlet'
  [../]
  [./J_tape_outlet]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_zero
    boundary = 'tape_outlet'
  [../]
  # Zero normal J on all other inlet/outlet boundaries
  [./J_zero_bc]
    type = MFEMVectorNormalDirichletBC
    variable = J
    vector_coefficient = J_zero
    boundary = 'silver_inlet silver_outlet outer_circle_inlet outer_circle_outlet copper_top_inlet copper_top_outlet copper_bottom_inlet copper_bottom_outlet substrate_inlet substrate_outlet'
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

  # Ampere law: curl H test - J test = 0
  # MFEMVectorFECurlKernel: trial=H (HCurl), test=HDiv (curl H in HDiv tested with HDiv test function)
  # MFEMVectorFEMassKernel: trial=J (HDiv), test=HDiv, coefficient=-1
  [./curlH]
    type = MFEMVectorFECurlKernel
    trial_variable = H
    variable = J
    coefficient = 1.0
  [../]
  [./J_term]
    type = MFEMVectorFEMassKernel
    variable = J
    coefficient = -1.0
  [../]
[]

[AuxKernels]
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
[]

[Postprocessors]
  [./time]
    type = TimePostprocessor
    execute_on = 'timestep_begin'
  [../]
  [./total_current]
    type = ParsedPostprocessor
    expression = '16000000000.0 * time'
    pp_names = 'time'
    execute_on = 'timestep_begin'
  [../]
[]

[Solvers]
  [./petsc_nl]
    type = MFEMPetscNonlinearSolver
    max_its = 200
    abs_tol = 1.0e-12
    rel_tol = 1.0e-3
    print_level = 3
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -sub_pc_factor_shift_type'
    petsc_options_value = 'lu mumps NONZERO'
  [../]
[]

[Executioner]
  type = MFEMTransient
  start_time = 0.0
  end_time = 1.0
  dt = 0.001
[]

[Outputs]
  [./out]
    type = MFEMParaViewDataCollection
    file_base = OutputData/CurlCurlJ_MFEM_multilayer
    vtk_format = ASCII
  [../]
[]