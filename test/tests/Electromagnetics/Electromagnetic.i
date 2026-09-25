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
  [./CurlComponent1]
    type = CurlComponent 
    H1 = Hz
    H2 = Hy
    component = 0
  [../]
  [./CurlComponent2]
    type = CurlComponent
    H1 = Hx
    H2 = Hz
    component = 1
  [../]
  [./CurlComponent3]
    type = CurlComponent
    H1 = Hy
    H2 = Hx
    component = 2
  [../]

  [./CurrentSource1]
    type = BodyForce
    variable = Jx
  [../]
  [./CurrentSource2]
    type = BodyForce
    variable = Jy
  [../]
  [./CurrentSource3]
    type = BodyForce
    variable = Jz
  [../]
[]

