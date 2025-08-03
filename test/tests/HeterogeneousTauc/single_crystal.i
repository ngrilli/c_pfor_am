# A modification of the CrystalPlasticityKalidindiUpdate
# to initialize tauc as spatially heterogeneous
# with a random value for each element.
# The purpose is to have a smoother yield point transition
# in the stress-strain curve.

[GlobalParams]
  displacements = 'ux uy uz'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 4
  elem_type = HEX8
[]

[AuxVariables]
  [fp_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [e_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [gss]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Physics/SolidMechanics/QuasiStatic/all]
  strain = FINITE
  add_variables = true
  generate_output = "strain_zz stress_zz"
[]

[AuxKernels]
  [fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = plastic_deformation_gradient
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
  [e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = total_lagrangian_strain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
  [gss]
    type = MaterialStdVectorAux
    variable = gss
    property = slip_resistance
    index = 0
    execute_on = timestep_end
  []
[]

[BCs]
  [symmy]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0
  []
  [symmx]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0
  []
  [symmz]
    type = DirichletBC
    variable = uz
    boundary = back
    value = 0
  []
  [tdisp]
    type = FunctionDirichletBC
    variable = uz
    boundary = front
    function = 0.001*t
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
  []
  [stress]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    maximum_substep_iteration = 10
  []
  [trial_xtalpl]
    type = CrystalPlasticityCopper
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
    resistance_tol = 1.0e-2
    # This parameter sets the random spatial variation of the initial gss
    # which is by default constant and uniform at 60.8 MPa
    gss_initial_std = 20.0 # MPa
    h = 1000.0
  []
[]

[Postprocessors]
  [stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  []
  [strain_zz]
    type = ElementAverageValue
    variable = strain_zz
  []
  [fp_zz]
    type = ElementAverageValue
    variable = fp_zz
  []
  [e_zz]
    type = ElementAverageValue
    variable = e_zz
  []
  [gss]
    type = ElementAverageValue
    variable = gss
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
  
  nl_abs_tol = 1e-8
  nl_rel_step_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_abs_step_tol = 1e-8
    
  dt = 0.01
  dtmax = 0.01
  dtmin = 0.001
  # run until time 50s to see stress-strain curve
  end_time = 0.01 #50.0
[]

[Outputs]
  exodus = true
  csv = true
[]
