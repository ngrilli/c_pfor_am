# Thermal sub-app for Tape geometry
# Solves for temperature on full domain using MFEMHypreBoomerAMG preconditioner
# Joule heating is received from J_field_app via MultiApp transfer
# Heat conductivity and specific heat are functions of temperature

[Mesh]
  type = MFEMFileMesh
  file = 'multilayer100um.msh'
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  # H1FESpace for temperature variable
  [./H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  [../]
  # L2FESpace for joule heating AuxVariable (received from J_field_app)
  [./L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = CONSTANT
  [../]
[]

[Variables]
  [./temperature]
    type = MFEMVariable
    fespace = H1FESpace
    initial_condition = 293.15  # Initial temperature 20°C
  [../]
[]

# Joule heating received from J_field_app via transfer
[AuxVariables]
  [./joule_heating]
    type = MFEMVariable
    fespace = L2FESpace
  [../]
  [./thermal_conductivity_out]
    type = MFEMVariable
    fespace = L2FESpace
  [../]
  [./specific_heat_out]
    type = MFEMVariable
    fespace = L2FESpace
  [../]
[]

[Functions]
  # Thermal conductivity as function of temperature (W/m·K)
  # Example: k(T) = k0 * (1 + alpha_k * (T - T0))
  [./thermal_conductivity_tape]
    type = MFEMParsedFunction
    expression = '20.0 * (1.0 + 0.003 * (T - 293.15))'  # Copper-like temperature dependence
    symbol_names = 'T'
    symbol_values = 'temperature'
  [../]
  [./thermal_conductivity_copper]
    type = MFEMParsedFunction
    expression = '400.0 * (1.0 - 0.001 * (T - 293.15))'  # Copper
    symbol_names = 'T'
    symbol_values = 'temperature'
  [../]
  [./thermal_conductivity_silver]
    type = MFEMParsedFunction
    expression = '429.0 * (1.0 - 0.001 * (T - 293.15))'  # Silver
    symbol_names = 'T'
    symbol_values = 'temperature'
  [../]
  [./thermal_conductivity_substrate]
    type = MFEMParsedFunction
    expression = '1.5'  # Substrate (constant)
  [../]
  [./thermal_conductivity_outer]
    type = MFEMParsedFunction
    expression = '0.5'  # Outer insulation (constant)
  [../]

  # Specific heat as function of temperature (J/kg·K)
  [./specific_heat_tape]
    type = MFEMParsedFunction
    expression = '385.0 * (1.0 + 0.0005 * (T - 293.15))'  # Copper-like
    symbol_names = 'T'
    symbol_values = 'temperature'
  [../]
  [./specific_heat_copper]
    type = MFEMParsedFunction
    expression = '385.0 * (1.0 + 0.0005 * (T - 293.15))'  # Copper
    symbol_names = 'T'
    symbol_values = 'temperature'
  [../]
  [./specific_heat_silver]
    type = MFEMParsedFunction
    expression = '235.0 * (1.0 + 0.0005 * (T - 293.15))'  # Silver
    symbol_names = 'T'
    symbol_values = 'temperature'
  [../]
  [./specific_heat_substrate]
    type = MFEMParsedFunction
    expression = '700.0'  # Substrate
  [../]
  [./specific_heat_outer]
    type = MFEMParsedFunction
    expression = '1000.0'  # Outer insulation
  [../]

  # Density (kg/m³) - constant for each material
  [./density_tape]
    type = MFEMParsedFunction
    expression = '8960.0'
  [../]
  [./density_copper]
    type = MFEMParsedFunction
    expression = '8960.0'
  [../]
  [./density_silver]
    type = MFEMParsedFunction
    expression = '10490.0'
  [../]
  [./density_substrate]
    type = MFEMParsedFunction
    expression = '2700.0'
  [../]
  [./density_outer]
    type = MFEMParsedFunction
    expression = '1200.0'
  [../]

  # Combined rho*Cp for time derivative
  [./rhoCp_tape]
    type = MFEMParsedFunction
    expression = 'density_tape * specific_heat_tape'
    symbol_names = 'density_tape specific_heat_tape'
    symbol_values = 'density_tape specific_heat_tape'
  [../]
  [./rhoCp_copper]
    type = MFEMParsedFunction
    expression = 'density_copper * specific_heat_copper'
    symbol_names = 'density_copper specific_heat_copper'
    symbol_values = 'density_copper specific_heat_copper'
  [../]
  [./rhoCp_silver]
    type = MFEMParsedFunction
    expression = 'density_silver * specific_heat_silver'
    symbol_names = 'density_silver specific_heat_silver'
    symbol_values = 'density_silver specific_heat_silver'
  [../]
  [./rhoCp_substrate]
    type = MFEMParsedFunction
    expression = 'density_substrate * specific_heat_substrate'
    symbol_names = 'density_substrate specific_heat_substrate'
    symbol_values = 'density_substrate specific_heat_substrate'
  [../]
  [./rhoCp_outer]
    type = MFEMParsedFunction
    expression = 'density_outer * specific_heat_outer'
    symbol_names = 'density_outer specific_heat_outer'
    symbol_values = 'density_outer specific_heat_outer'
  [../]
[]

[FunctorMaterials]
  # Thermal conductivity materials
  [./thermal_conductivity_tape_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'thermal_conductivity'
    prop_values = 'thermal_conductivity_tape'
    block = 'tape'
  [../]
  [./thermal_conductivity_copper_bottom_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'thermal_conductivity'
    prop_values = 'thermal_conductivity_copper'
    block = 'copper_bottom'
  [../]
  [./thermal_conductivity_copper_top_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'thermal_conductivity'
    prop_values = 'thermal_conductivity_copper'
    block = 'copper_top'
  [../]
  [./thermal_conductivity_silver_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'thermal_conductivity'
    prop_values = 'thermal_conductivity_silver'
    block = 'silver'
  [../]
  [./thermal_conductivity_substrate_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'thermal_conductivity'
    prop_values = 'thermal_conductivity_substrate'
    block = 'substrate'
  [../]
  [./thermal_conductivity_outer_circle_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'thermal_conductivity'
    prop_values = 'thermal_conductivity_outer'
    block = 'outer_circle'
  [../]

  # rho*Cp materials for time derivative
  [./rhoCp_tape_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'rhoCp'
    prop_values = 'rhoCp_tape'
    block = 'tape'
  [../]
  [./rhoCp_copper_bottom_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'rhoCp'
    prop_values = 'rhoCp_copper'
    block = 'copper_bottom'
  [../]
  [./rhoCp_copper_top_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'rhoCp'
    prop_values = 'rhoCp_copper'
    block = 'copper_top'
  [../]
  [./rhoCp_silver_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'rhoCp'
    prop_values = 'rhoCp_silver'
    block = 'silver'
  [../]
  [./rhoCp_substrate_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'rhoCp'
    prop_values = 'rhoCp_substrate'
    block = 'substrate'
  [../]
  [./rhoCp_outer_circle_mat]
    type = MFEMGenericFunctorMaterial
    prop_names = 'rhoCp'
    prop_values = 'rhoCp_outer'
    block = 'outer_circle'
  [../]
[]

[Kernels]
  # Heat conduction: k * grad(T) · grad(test)
  [./diffusion]
    type = MFEMDiffusionKernel
    variable = temperature
    coefficient = thermal_conductivity
  [../]

  # Transient term: rho*Cp * dT/dt
  [./dT_dt]
    type = MFEMTimeDerivativeMassKernel
    variable = temperature
    coefficient = rhoCp
  [../]

  # Joule heating source term
  [./joule_source]
    type = MFEMDomainLFKernel
    variable = temperature
    coefficient = joule_heating
  [../]
[]

[BCs]
  # Convective cooling on outer boundary
  [./outer_convection]
    type = MFEMConvectiveHeatFluxBC
    variable = temperature
    boundary = 'outer_circle_external'
    T_infinity = 293.15
    heat_transfer_coefficient = 10.0
  [../]

  # Symmetry/adiabatic on other external boundaries (implicit natural BC)
[]

[AuxKernels]
  # Output thermal conductivity
  [./thermal_conductivity_aux]
    type = MFEMScalarProjectionAux
    variable = thermal_conductivity_out
    coefficient = thermal_conductivity
  [../]

  # Output specific heat
  [./specific_heat_aux]
    type = MFEMScalarProjectionAux
    variable = specific_heat_out
    coefficient = specific_heat_tape  # Note: this will show tape value everywhere; use block-restricted for per-block
  [../]
[]

[Solvers]
  [./boomeramg]
    type = MFEMHypreBoomerAMG
  [../]
  [./main]
    type = MFEMHypreGMRES
    preconditioner = boomeramg
    l_tol = 1e-10
    l_max_its = 1000
    print_level = 1
  [../]
  [./newton]
    type = MFEMNewtonNonlinearSolver
    use_initial_guess = true
    max_its = 20
    abs_tol = 1e-8
    rel_tol = 1e-6
    print_level = 1
  [../]
[]

[Executioner]
  type = MFEMTransient
  start_time = 0.0
  end_time = 1.11
  dt = 0.0025
[]

[Outputs]
  [./out]
    type = MFEMParaViewDataCollection
    file_base = OutputData/tape_thermal_submesh
    vtk_format = ASCII
    time_step_interval = 40
  [../]
[]