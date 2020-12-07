[Mesh]
  [./gen_mesh]
    type = FileMeshGenerator
    file = DummyVar_out.e
  [../]
  [./add_block_id_to_delete]
    type = ElementSubdomainIDGenerator
	input = gen_mesh
    element_ids = '6 7 8 11 12 13 16 17 18'
    subdomain_ids = '2 2 2 2 2 2 2 2 2'
  [../]
  [./deleteelement]
    type = BlockDeletionGenerator
    block_id = '2'
    input = add_block_id_to_delete
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = disp_x_ic
    [../]
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = disp_y_ic
    [../]
  [../]

  [./disp_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = disp_z_ic
    [../]
  [../]
  
  [./dummy_var]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = dummy_var_ic
    [../]
  [../]  
[]

[AuxVariables]

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./dummy_aux_var]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = dummy_aux_var_ic
    [../]
  [../]

[]

[Functions]
  [./disp_x_ic]
    type = SolutionFunction
    solution = soln_disp_x
  [../]
  [./disp_y_ic]
    type = SolutionFunction
    solution = soln_disp_y
  [../]
  [./disp_z_ic]
    type = SolutionFunction
    solution = soln_disp_z
  [../]
  [./dummy_var_ic]
    type = SolutionFunction
    solution = soln_dummy_var
  [../]
  [./dummy_aux_var_ic]
    type = SolutionFunction
    solution = soln_dummy_aux_var
  [../]
[]

[UserObjects]
  [./soln_disp_x]
    type = SolutionUserObject
    mesh = DummyVar_out.e
	timestep = 'LATEST'
    system_variables = disp_x
  [../]
  [./soln_disp_y]
    type = SolutionUserObject
    mesh = DummyVar_out.e
	timestep = 'LATEST'
    system_variables = disp_y
  [../]
  [./soln_disp_z]
    type = SolutionUserObject
    mesh = DummyVar_out.e
	timestep = 'LATEST'
    system_variables = disp_z
  [../]
  [./soln_dummy_var]
    type = SolutionUserObject
    mesh = DummyVar_out.e
	timestep = 'LATEST'
    system_variables = dummy_var
  [../]
  [./soln_dummy_aux_var]
    type = SolutionUserObject
    mesh = DummyVar_out.e
	timestep = 'LATEST'
    system_variables = dummy_aux_var
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    add_variables = true
  [../]
  [./ddummydt]
    type = TimeDerivative
    variable = dummy_var
  [../]
[]

[AuxKernels]

  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]  

  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]

  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]

[]

[BCs]
  [./z0]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./y0]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./x0]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  [./x1]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0.0
  [../]
[]
 
[Postprocessors]

[]

[Materials]
  [./stress]
    type = ComputeFiniteStrainElasticStress
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '1.5e6 0.75e6 0.75e6 1.5e6 0.75e6 1.5e6 0.375e6 0.375e6 0.375e6'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
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
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg          31'
  line_search = 'none'
  l_max_its = 50
  nl_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  l_tol = 1e-8

  start_time = 0.01
  end_time = 0.02
  dt = 0.005
  dtmin = 0.00001
[]

[Outputs]
  csv = false
  [./my]
    type = Checkpoint
    num_files = 2
    interval = 1
  [../]
  [./out]
    type = Exodus
	interval = 1
	checkpoint = true
  [../]
[]
