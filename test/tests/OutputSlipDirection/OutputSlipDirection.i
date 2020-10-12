[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  elem_type = HEX8
  displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    # Use the initial Condition block underneath the variable
    # for which we want to apply this initial condition
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./temp]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = '293.0'
    [../]
  [../]

  [./slipdirection1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slipdirection2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slipdirection3]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[Functions]
  [./temperature_load]
    type = ParsedFunction
    value = '293.0'
  [../]

  [./disp_load]
    type = ParsedFunction
    value = '0.1*t'
  [../]
[]

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'euler_ang_test.inp'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    ngrain = 1
    read_type = indexgrain
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    add_variables = true
  [../]
[]

[AuxKernels]

# check how the slip direction
# 1	0  -1
# rotates under a rotation with Euler angles
# -45.0 0.0 0.0
  [./slipdirection1]
    type = MaterialStdVectorAux
    variable = slipdirection1
    property = slip_direction
    index = 21
    execute_on = timestep_end
  [../]
  [./slipdirection2]
    type = MaterialStdVectorAux
    variable = slipdirection2
    property = slip_direction
    index = 22
    execute_on = timestep_end
  [../]
  [./slipdirection3]
    type = MaterialStdVectorAux
    variable = slipdirection3
    property = slip_direction
    index = 23
    execute_on = timestep_end
  [../]

[]

[BCs]
  [./z0_bot]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./y0_bot]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./x0_bot]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  [./z1_disp]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = disp_load
  [../]
[]
 
[Postprocessors]

[]


[Materials]
  [./crysp]
    type = FiniteStrainCrystalPlasticityThermal
    gtol = 1e-2
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 12 #Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' # slip rate equations parameters 
    hprops = '1.0 541.5 136.0 200.0 2.5' # hardening properties
    gprops = '1 12 136.0' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
    thermal_expansion = '17e-6'
    reference_temperature = '293.0'
    temp = temp
# Calibrated using table 1 in:
# M.R. DAYMOND and P.J. BOUCHARD
# Elastoplastic Deformation of 316 Stainless Steel Under
# Tensile Loading at Elevated Temperatures
# METALLURGICAL AND MATERIALS TRANSACTIONS A
# VOLUME 37A, JUNE 2006—1873
	dCRSS_dT_A = 0.53
	dCRSS_dT_B = 0.47
	dCRSS_dT_C = 0.008
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPGrain
# Elastic constants of 316L SS from:
# Clausen, B., Lorentzen, T. and Leffers, T. 
# Self-consistent modelling of the plastic
# deformation of FCC polycrystals and its implications for diffraction
# measurements of internal stresses.
# Acta Mater. 46, 3087–3098 (1998).
    C_ijkl = '2.046e5 1.377e5 1.377e5 2.046e5 1.377e5 2.046e5 1.262e5 1.262e5 1.262e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
    temp = temp
    dC11_dT = 0.0004757
    dC12_dT = 0.0004757
    dC44_dT = 0.0004757
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

  start_time = 0.0
  end_time = 0.02
  dt = 0.01
  dtmin = 0.01
[]

[Outputs]
  csv = true
  [./out]
    type = Exodus
  [../]
[]
