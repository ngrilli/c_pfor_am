[Mesh]
  [abaqus_mesh]
    type = FileMeshGenerator
    file = 'Job-0_0_0_121.inp'
  []
[]

[Variables]
  [./u_x]
  [../]
  [./u_y]
  [../]
  [./u_z]
  [../]
  [./global_strain]
    order = SIXTH
    family = SCALAR
  [../]
[]

[AuxVariables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./s00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e22]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./disp_x]
    type = GlobalDisplacementAux
    variable = disp_x
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 0
  [../]
  [./disp_y]
    type = GlobalDisplacementAux
    variable = disp_y
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 1
  [../]
  [./disp_z]
    type = GlobalDisplacementAux
    variable = disp_z
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 2
  [../]
  [./s00]
    type = RankTwoAux
    variable = s00
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  [../]
  [./s11]
    type = RankTwoAux
    variable = s11
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  [../]  
  [./s22]
    type = RankTwoAux
    variable = s22
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  [../]
  [./e00]
    type = RankTwoAux
    variable = e00
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 0
  [../]
  [./e11]
    type = RankTwoAux
    variable = e11
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 1
  [../]
  [./e22]
    type = RankTwoAux
    variable = e22
    rank_two_tensor = total_strain
    index_i = 2
    index_j = 2
  [../]
[]

[GlobalParams]
  displacements = 'u_x u_y u_z'
[]

[Kernels]
  [./TensorMechanics]
  [../]
[]

[ScalarKernels]
  [./global_strain]
    type = GlobalStrain
    variable = global_strain
    global_strain_uo = global_strain_uo
  [../]
[]

[BCs]

  [./Periodic]
  
    [./x]
      variable = u_x
      primary = face_x_minus
      secondary = face_x_plus
      translation = '1 0 0'
    [../]  
    
    [./y]
      variable = u_y
      primary = face_y_minus
      secondary = face_y_plus
      translation = '0 1 0'
    [../]  
    
    [./z]
      variable = u_z
      primary = face_z_minus
      secondary = face_z_plus
      translation = '0 0 1'
    [../]  
 
  [../]

  # fix center point location
  [./centerfix_x]
    type = DirichletBC
    boundary = origin
    variable = u_x
    value = 0
  [../]
  [./centerfix_y]
    type = DirichletBC
    boundary = origin
    variable = u_y
    value = 0
  [../]
  [./centerfix_z]
    type = DirichletBC
    boundary = origin
    variable = u_z
    value = 0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '70e9 0.33'
    fill_method = symmetric_isotropic_E_nu
  [../]
  [./strain]
    type = ComputeSmallStrain
    global_strain = global_strain
  [../]
  [./global_strain]
    type = ComputeGlobalStrain
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[UserObjects]
  [./global_strain_uo]
    type = GlobalStrainUserObject
    applied_stress_tensor = '0 0 0 5e9 5e9 5e9'
    execute_on = 'Initial Linear Nonlinear'
  [../]
[]

[Postprocessors]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'

  line_search = basic

  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm         31   preonly   lu      1'

  l_max_its = 30
  nl_max_its = 12

  nl_rel_tol = 1.0e-10

  start_time = 0.0
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
