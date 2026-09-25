[Mesh]
  [circle]
    type = ConcentricCircleMeshGenerator
    has_outer_square = false
    radii = 1
    num_sectors = 10
    rings = 1
    preserve_volumes = false
  []
  [side]
    type = SideSetsAroundSubdomainGenerator
    input = circle
    new_boundary = side
    block = 1
  []
  [./extrude]
    input = side
    type = MeshExtruderGenerator
    num_layers = 20
    extrusion_vector = '0 0 2'
    bottom_sideset = 'bottom'
    top_sideset = 'top'
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

  # Faraday law
  [./CurlEComponent1]
    type = PowerLawCurlComponent
    variable = Jx
    J1 = Jz
    J2 = Jy
    E0 = 1.0
    J0 = 1.0
    n = 21.0
    component = 0
  [../]
  [./CurlEComponent2]
    type = PowerLawCurlComponent
    variable = Jy
    J1 = Jx
    J2 = Jz
    E0 = 1.0
    J0 = 1.0
    n = 21.0
    component = 1
  [../]
  [./CurlEComponent3]
    type = PowerLawCurlComponent
    variable = Jz
    J1 = Jy
    J2 = Jx
    E0 = 1.0
    J0 = 1.0
    n = 21.0
    component = 2
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
    boundary = side
    value = 0.0
  [../]
  [./sideHy]
    type = DirichletBC
    variable = Hy
    boundary = side
    value = 0.0
  [../]
  [./sideHz]
    type = DirichletBC
    variable = Hz
    boundary = side
    value = 0.0
  [../]
  [./topJz]
    type = DirichletBC
    variable = Jz
    boundary = top
    value = 1.0
  [../]
  [./bottomJz]
    type = DirichletBC
    variable = Jz
    boundary = bottom
    value = 1.0
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
  
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-5
  
  start_time = 0.0
  end_time = 0.1
  dt = 0.001

[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]
