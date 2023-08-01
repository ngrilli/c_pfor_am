# SubApp input file for phase field fracture
# staggered solver

# same mesh as master app

[Mesh]
  [./nepermesh]
    type = FileMeshGenerator
    file = 'n8-id1.msh'
  [../]
  [./left_modifier]
    type = BoundingBoxNodeSetGenerator
    input = nepermesh
    new_boundary = left
    top_right = '0.001 40.001 1.001'
    bottom_left = '-0.001 -0.001 -0.001'
  [../]
  [./bottom_modifier]
    type = BoundingBoxNodeSetGenerator
    input = left_modifier
    new_boundary = bottom
    top_right = '40.001 0.001 1.001'
    bottom_left = '-0.001 -0.001 -0.001'
  [../]
  [./back_modifier]
    type = BoundingBoxNodeSetGenerator
    input = bottom_modifier
    new_boundary = back
    top_right = '40.001 40.001 0.001'
    bottom_left = '-0.001 -0.001 -0.001'
  [../]
  [./right_modifier]
    type = BoundingBoxNodeSetGenerator
    input = back_modifier
    new_boundary = right
    top_right = '40.001 40.001 1.001'
    bottom_left = '39.999 -0.001 -0.001'
  [../]
[]

[AuxVariables]
  # set bound variable to impose
  # that damage can only grow
  [./bounds_dummy]
  [../]
  
  [./elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./delastic_energy_dc]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./d2elastic_energy_dc2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
[]

[Modules/PhaseField/Nonconserved/c]
  free_energy = F
  kappa = kappa_op
  mobility = L
[]

# bounds used to impose that damage can only grow
[Bounds]
  [./irreversibility]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = lower
  [../]
  [./upper]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = upper
    bound_value = 1.0
  [../]
[]


[Materials]

  # gc_prop is the fracture energy
  # if plastic work is included
  # fracture takes place at lower strain
  # l is the characteristic crack diffusion length
  # ideally element size should be about one fourth of l
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '50.0 1.0 1e-3'
  [../]

  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0/(gc_prop * visco)'
  [../]

  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l'
  [../]
  
  [./elastic_energy]
    type = DerivativeParsedMaterial
    f_name = elastic_energy
    args = 'elastic_energy'
    function = 'elastic_energy'
    derivative_order = 2
  [../]
  
  [./delastic_energy_dc]
    type = DerivativeParsedMaterial
    f_name = delastic_energy/dc
    args = 'delastic_energy_dc'
    function = 'delastic_energy_dc'
    derivative_order = 1
  [../]
  
  [./d2elastic_energy_dc2]
    type = DerivativeParsedMaterial
    f_name = d^2elastic_energy/dc^2
    args = 'd2elastic_energy_dc2'
    function = 'd2elastic_energy_dc2'
    derivative_order = 0
  [../]

  [./local_fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'c'
    material_property_names = 'gc_prop l'
    function = 'c^2 * gc_prop / 2 / l'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    args = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    f_name = F
  [../]

[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  # variational inequality solver is needed to avoid damage to grow larger than 1
  petsc_options_iname = '-pc_type  -snes_type'
  petsc_options_value = 'lu vinewtonrsls'

  line_search = 'none'
  
  nl_abs_tol = 1e-3
  nl_rel_tol = 1e-3

  nl_max_its = 10
  l_max_its = 50
  
  [./TimeStepper]
    type = ConstantDT
    dt = 0.01
    growth_factor = 10
  [../]

  dtmax = 0.01
  dtmin = 1.0e-40
  end_time = 0.01 # run until 100 s to see full crack propagation
[]

[Outputs]
  [./out]
    type = Exodus
    interval = 1 #100
  [../]
  # uncomment for restart option
  #[./restart]
  #  type = Checkpoint
  #  num_files = 5
  #  interval = 10
  #[../]
[]

