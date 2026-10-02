[Mesh]
  [./generated_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 20
    ny = 20
    nz = 20
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 1.0
  [../]
  [./add_wire]
    type = ParsedSubdomainMeshGenerator
    input = generated_mesh
    combinatorial_geometry = '(x-0.5)^2 < 0.1*0.1 & (y-0.5)^2 < 0.1*0.1'
    block_id = 1
    block_name = 'wire'
  [../]
  [./current_density_inlet]
    type = BoundingBoxNodeSetGenerator
    input = add_wire
    bottom_left = '0.4 0.4 0.0'
    top_right = '0.6 0.6 0.0'
    new_boundary = 'current_density_inlet'
  [../]
[]

[Variables]
  [./Hx]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Hy]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Hz]
    order = FIRST
    family = LAGRANGE
  [../]

  [./Jx]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Jy]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Jz]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]

  # Ampere law
  [./CurlHComponent1]
    type = CurlComponent 
    variable = Hx
    H1 = Hz
    H2 = Hy
    component = 0
  [../]
  [./CurlHComponent2]
    type = CurlComponent
    variable = Hy
    H1 = Hx
    H2 = Hz
    component = 1
  [../]
  [./CurlHComponent3]
    type = CurlComponent
    variable = Hz
    H1 = Hy
    H2 = Hx
    component = 2
  [../]

  [./CurrentSource1]
    type = CoupledForce
    variable = Hx
    v = Jx
  [../]
  [./CurrentSource2]
    type = CoupledForce
    variable = Hy
    v = Jy
  [../]
  [./CurrentSource3]
    type = CoupledForce
    variable = Hz
    v = Jz
  [../]

  # Faraday law air
  [./CurlEComponent1_air]
    type = PowerLawCurlComponent
    variable = Jx
    J1 = Jz
    J2 = Jy
    E0 = 1.0
    J0 = 1.0
    n = 1.0
    component = 0
    block = 0
  [../]
  [./CurlEComponent2_air]
    type = PowerLawCurlComponent
    variable = Jy
    J1 = Jx
    J2 = Jz
    E0 = 1.0
    J0 = 1.0
    n = 1.0
    component = 1
    block = 0
  [../]
  [./CurlEComponent3_air]
    type = PowerLawCurlComponent
    variable = Jz
    J1 = Jy
    J2 = Jx
    E0 = 1.0
    J0 = 1.0
    n = 1.0
    component = 2
    block = 0
  [../]

  # Faraday law superconductor
  [./CurlEComponent1_superconductor]
    type = PowerLawCurlComponent
    variable = Jx
    J1 = Jz
    J2 = Jy
    E0 = 1.0
    J0 = 1.0
    n = 21.0
    component = 0
    block = 1
  [../]
  [./CurlEComponent2_superconductor]
    type = PowerLawCurlComponent
    variable = Jy
    J1 = Jx
    J2 = Jz
    E0 = 1.0
    J0 = 1.0
    n = 21.0
    component = 1
    block = 1
  [../]
  [./CurlEComponent3_superconductor]
    type = PowerLawCurlComponent
    variable = Jz
    J1 = Jy
    J2 = Jx
    E0 = 1.0
    J0 = 1.0
    n = 21.0
    component = 2
    block = 1
  [../]

  [./TimeDerivative1]
    type = CoupledTimeDerivative
    variable = Jx
    v = Hx
  [../]
  [./TimeDerivative2]
    type = CoupledTimeDerivative
    variable = Jy
    v = Hy
  [../]
  [./TimeDerivative3]
    type = CoupledTimeDerivative
    variable = Jz
    v = Hz
  [../]
[]

[BCs]
  [./sideHx]
    type = DirichletBC
    variable = Hx
    boundary = 'left right top bottom'
    value = 0.0
  [../]
  [./sideHy]
    type = DirichletBC
    variable = Hy
    boundary = 'left right top bottom'
    value = 0.0
  [../]
  [./sideHz]
    type = DirichletBC
    variable = Hz
    boundary = 'left right top bottom'
    value = 0.0
  [../]
  [./current_density_inlet_BC]
    type = FunctionDirichletBC
    variable = Jz
    boundary = current_density_inlet
    function = ramp_up_current_density
  [../]
[]

[Functions]
  [./ramp_up_current_density]
    type = PiecewiseLinear
    x = '0.0 0.1 1.0' # time
    y = '0.0 0.1 0.1' # current density
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
  solve_type = 'NEWTON'
  
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg          51'
  #line_search = 'none'
  
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  
  start_time = 0.0
  end_time = 0.1
  dt = 0.001

[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]
