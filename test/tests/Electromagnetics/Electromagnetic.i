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