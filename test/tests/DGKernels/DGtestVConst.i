[Mesh]

  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 1
    xmax = 2.0
    ymax = 2.0
    zmax = 0.2
    elem_type = HEX8 
  [../] 

[]

[GlobalParams]

[]

[Variables]
 
  [./rho_t]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = init_rho_t
    [../]
  [../]
  
  [./rho_gnd_edge]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = init_rho_gnd_edge
    [../]
  [../]
  
[]

[AuxVariables]

[]

[Functions]

  [./dts]
    type = PiecewiseConstant
    x = '0.0 10.0'
    y = '0.1 0.1'
  [../]
  
  [./init_rho_gnd_edge]
    type = ParsedFunction
	value = 'if(3.0*(x-1.0),0.0,1.0)'
  [../]
  
  [./init_rho_t]
    type = ParsedFunction
	value = 'if(3.0*(x-1.0),0.0,1.0)'
  [../]
  
[]

[UserObjects]

[]

[Kernels]

  [./drho_gnd_edge_dt]
    type = TimeDerivative
    variable = rho_gnd_edge
  [../]

  [./drho_t_dt]
    type = TimeDerivative
    variable = rho_t
  [../]  
    
[]

[DGKernels]

  [./rho_t_advection_edge]
    type = DGAdvectionCoupledVConst
    variable = rho_gnd_edge
	rho_coupled = rho_t
    velocity = '1 0 0'
  [../] 
  
  [./rho_gnd_edge_advection]
    type = DGAdvectionCoupledVConst
    variable = rho_t
	rho_coupled = rho_gnd_edge
    velocity = '1 0 0'
  [../]

[]

[AuxKernels]

[]

[BCs]
  
  [./Periodic]
  
    [./auto_rho_t_boundary_x]
      variable = rho_t
      primary = 'left'
	  secondary = 'right'
	  translation = '2.0 0.0 0.0'
    [../]
    [./auto_rho_gnd_edge_boundary_x]
      variable = rho_gnd_edge
      primary = 'left'
	  secondary = 'right'
      translation = '2.0 0.0 0.0'  
    [../]
	
    [./auto_rho_t_boundary_y]
      variable = rho_t
      primary = 'bottom'
	  secondary = 'top'
	  translation = '0.0 2.0 0.0'
    [../]
    [./auto_rho_gnd_edge_boundary_y]
      variable = rho_gnd_edge
      primary = 'bottom'
	  secondary = 'top'
      translation = '0.0 2.0 0.0'	  
    [../]
	
    [./auto_rho_t_boundary_z]
      variable = rho_t
      primary = 'back'
	  secondary = 'front'
	  translation = '0.0 0.0 0.2'
    [../]
    [./auto_rho_gnd_edge_boundary_z]
      variable = rho_gnd_edge
      primary = 'back'
	  secondary = 'front'
	  translation = '0.0 0.0 0.2'
    [../]

  [../]  
  
[]
 
[Postprocessors]

[]

[Materials]

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
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  
  line_search = 'none'
  automatic_scaling = true
  
  l_max_its = 50
  nl_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  
  [./TimeStepper]
    type = FunctionDTGrowth
    function = dts
	cutback_factor_at_failure = 0.1
	growth_factor = 1.2
  [../]

  start_time = 0.0
  end_time = 0.2 #10.0
  
  dtmin = 1.0e-10
  timestep_tolerance = 1.0e-10
[]

[Outputs]

  csv = false

  [./out]
    type = Exodus
    interval = 1
  [../]
  
[]
