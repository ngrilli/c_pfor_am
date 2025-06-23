[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Mesh]
  [./file_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 20
    ny = 20
    nz = 8
    xmin = -0.00020
    xmax = 0.00020
    ymin = -0.00020
    ymax = 0.00020
    zmin = 0.00000
    zmax = 0.00016
  [../]
  
  # separate subdomain that include printed layers and inactive subdomain
  [./top_layer]
    type = SubdomainBoundingBoxGenerator
    input = file_mesh
    bottom_left = '-0.00020 -0.00020 0.00012'
    top_right = '0.00020 0.00020 0.00016'
    block_id = 1
  [../]
  
  # separate subdomains for first and second printed layers
  [./first_layer_center]
    type = SubdomainBoundingBoxGenerator
    input = top_layer
    bottom_left = '-0.00010 -0.00010 0.00012'
    top_right = '0.00010 0.00010 0.00014'
    block_id = 2
  [../]
  [./second_layer_center]
    type = SubdomainBoundingBoxGenerator
    input = first_layer_center
    bottom_left = '-0.00010 -0.00010 0.00014'
    top_right = '0.00010 0.00010 0.00016'
    block_id = 3
  [../]  
  
  # define sidesets on substrate for boundary conditions
  [./front_face_sideset]
    input = second_layer_center
    type = SideSetsAroundSubdomainGenerator
    normal = '0 0 1'
    block = 0
    new_boundary = 'substrate_front'
  [../]
  [./left_face_sideset]
    input = front_face_sideset
    type = SideSetsAroundSubdomainGenerator
    normal = '-1 0 0'
    block = 0
    new_boundary = 'substrate_left'
  [../]
  [./right_face_sideset]
    input = left_face_sideset
    type = SideSetsAroundSubdomainGenerator
    normal = '1 0 0'
    block = 0
    new_boundary = 'substrate_right'
  [../]  
  [./bottom_face_sideset]
    input = right_face_sideset
    type = SideSetsAroundSubdomainGenerator
    normal = '0 -1 0'
    block = 0
    new_boundary = 'substrate_bottom'
  [../]  
  [./top_face_sideset]
    input = bottom_face_sideset
    type = SideSetsAroundSubdomainGenerator
    normal = '0 1 0'
    block = 0
    new_boundary = 'substrate_top'
  [../]    
  
[]

# Every parameter in the GlobalParams block will be inserted into every block/sub-block where that parameter name is defined
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

# temperature variable
[Variables]
  [./temp]
    block = 0
    initial_condition = 293.0
  [../]
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./all]
        strain = SMALL
        incremental = true
        add_variables = true
        eigenstrain_names = eigenstrain
        generate_output = 'mechanical_strain_xx mechanical_strain_xy mechanical_strain_xz mechanical_strain_yy mechanical_strain_yz mechanical_strain_zz elastic_strain_xx elastic_strain_xy elastic_strain_xz elastic_strain_yy elastic_strain_yz elastic_strain_zz plastic_strain_xx plastic_strain_xy plastic_strain_xz plastic_strain_yy plastic_strain_yz plastic_strain_zz stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz'
        use_automatic_differentiation = false
        temperature = temp
        volumetric_locking_correction = true
        block = 0
      [../]
    [../]
  [../]
[]

# heat equation
[Kernels]
  [./time]
    type = HeatConductionTimeDerivative
    variable = temp
    block = 0
  [../]
  [./heat_conduct]
    type = HeatConduction
    variable = temp
    block = 0
  [../]
  [./heat_source]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    block = 0
  [../]
[]

[AuxVariables]
  # auxiliary variable to determine element activation
  [level_set]
  []
  
  # auxiliary variable to visualize heat source 
  [volumetric_heat]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  []
[]

[AuxKernels]
  [level_set]
    type = FunctionAux
    variable = level_set
    function = level_set
  []
  [volumetric_heat]
    type = ADMaterialRealAux
    variable = volumetric_heat
    property = volumetric_heat
    block = 0
  []
[]

[BCs]
  # Convective heat loss on front surface of the substrate
  # assuming other surfaces are connected to a larger substrate
  [./convective_front]
    type = ConvectiveHeatFluxBC
    variable = temp
    boundary = 'substrate_front'
    T_infinity = 293.0
    heat_transfer_coefficient = 20.0
  [../]

  # Radiative heat loss on front surface of the substrate
  # assuming other surfaces are connected to a larger substrate
  [./radiative_front]
    type = FunctionRadiativeBC
    variable = temp
    boundary = 'substrate_front'
    emissivity_function = '0.24'
    Tinfinity = 293.0
  [../]
  
  # The larger substrate is modeled as a heat sink
  [./heat_sink]
    type = DirichletBC
    variable = temp
    boundary = 'substrate_left substrate_right substrate_top substrate_bottom back'
    value = 293.0
  [../]

  # Mechanical boundary conditions on the substrate
  [./left_block_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'substrate_left substrate_right'
    value = 0.0
  [../]
  [./bottom_block_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'substrate_bottom substrate_top'
    value = 0.0
  [../]
  [./back_block_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  [../]
[]

[Materials]

  [./heat]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = k_function
    specific_heat_temperature_function = cp_function
    temp = temp
    block = 0
  [../]

  # Ti6Al4V density
  [./density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '4430'
    block = 0
  [../]

  # Laser parameters, see Table 1 in
  # Cook, Ritchie
  # Determining the laser absorptivity of Ti-6Al-4V during laser powder bed fusion by calibrated melt pool simulation
  # Optics & Laser Technology, Volume 162, July 2023, 109247
  # https://www.sciencedirect.com/science/article/pii/S0030399223001408
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    rx = 60e-6
    ry = 60e-6
    rz = 50e-6
    power = 100
    efficiency = 0.2
    factor = 1
    function_x = path_x
    function_y = path_y
    function_z = path_z
    block = 0
  [../]

  [./stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plas'
    max_iterations = 10000
    relative_tolerance = 1e-5
    absolute_tolerance = 1e-5
    combined_inelastic_strain_weights = '1.0'
    internal_solve_full_iteration_history = true
    perform_finite_strain_rotations = false
    tangent_operator = nonlinear
    block = 0
  [../]

  [./plas]
    type = IsotropicPlasticityStressUpdate
    yield_stress_function = yield_stress_function
    hardening_function = hardening_function
    temperature = temp
    block = 0
  [../]

  [./elasticity_tensor]
    type = ComputeVariableIsotropicElasticityTensor
    args = temp
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
    block = 0
  [../]

  # Coefficient of thermal expansion as a function of temperature
  # implementation based on modules/solid_mechanics/test/tests/thermal_expansion_function/finite_const.i
  [./thermal_expansion_strain]
    type = ComputeInstantaneousThermalExpansionFunctionEigenstrain
    thermal_expansion_function = alpha_function
    stress_free_temperature = 293 
    temperature = temp
    eigenstrain_name = eigenstrain
    block = 0
  [../]

  # Young modulus as a function of temperature
  # implementation based on modules/combined/test/tests/thermo_mech/youngs_modulus_function_temp.i
  # Table 2 in
  # Jun Cao, Michael A. Gharghouri, Philip Nash
  # Finite-element analysis and experimental validation of thermal residual stress and distortion in electron beam additive manufactured Ti-6Al-4V build plates
  # Journal of Materials Processing Technology, Volume 237, November 2016, Pages 409-419
  # https://www.sciencedirect.com/science/article/pii/S0924013616302126
  [./youngs_modulus]
    type = PiecewiseLinearInterpolationMaterial
    x = '293   400    500   600   700   800   900   1000  1100  1200  1300  1400  1950 19500'
    y = '102e9 101e9  95e9  91e9  85e9  80e9  75e9  70e9  65e9  60e9  35e9  20e9  10e9 10e9'
    property = youngs_modulus
    variable = temp
    block = 0
  [../]

  # Poisson ratio as a function of temperature
  # Table 2 in
  # Jun Cao, Michael A. Gharghouri, Philip Nash
  # Finite-element analysis and experimental validation of thermal residual stress and distortion in electron beam additive manufactured Ti-6Al-4V build plates
  # Journal of Materials Processing Technology, Volume 237, November 2016, Pages 409-419
  # https://www.sciencedirect.com/science/article/pii/S0924013616302126
  [./poissons_ratio]
    type = PiecewiseLinearInterpolationMaterial
    x = '293    400     500    600    700    800    900   1000   1100    1200   1300  1400    1950'
    y = '0.345 0.350 0.355 0.360 0.365 0.370 0.375 0.385 0.395 0.405 0.430 0.430 0.430'
    property = poissons_ratio
    variable = temp
    block = 0
  [../]

[]

[Functions]

  # path files: x, y, z coordinates of the heat source as a function of time
  [./path_x]
    type = PiecewiseLinear
    data_file = path_x.txt
    format = columns
  [../]
  [./path_y]
    type = PiecewiseLinear
    data_file = path_y.txt
    format = columns
  [../]
  [./path_z]
    type = PiecewiseLinear
    data_file = path_z.txt
    format = columns
  [../]

  # temperature dependent thermal conductivity
  # Jun Cao, Michael A. Gharghouri, Philip Nash
  # Finite-element analysis and experimental validation of thermal residual stress and distortion in electron beam additive manufactured Ti-6Al-4V build plates
  # Journal of Materials Processing Technology, Volume 237, November 2016, Pages 409-419
  # https://www.sciencedirect.com/science/article/pii/S0924013616302126
  [./k_function]
    type = PiecewiseLinear
    x = '298  373   473   573    673   773  873 973  1073 1173 1268 1268.01 1373 1473 1573 1673 1773 1873 1923 1923.01 1973'
    y = '7    7.45  8.75  10.15 11.35 12.6 14.2 15.5 17.8 20.2 22.7 19.3 21   22.9 23.7 24.6 25.8 27   28.4 33.4 34.6'
  [../]

  # temperature dependent specific heat
  # Jun Cao, Michael A. Gharghouri, Philip Nash
  # Finite-element analysis and experimental validation of thermal residual stress and distortion in electron beam additive manufactured Ti-6Al-4V build plates
  # Journal of Materials Processing Technology, Volume 237, November 2016, Pages 409-419
  # https://www.sciencedirect.com/science/article/pii/S0924013616302126
  [./cp_function]
    type = PiecewiseLinear
    x = '298 373 473 573 673 773 873 973 1073 1173 1268 1268.01 1373 1473 1573 1673 1773 1873 1923 1923.01 1973 2073 2173'
    y = '546 562 584 606 629 651 673 694 714  734  753  641  660  678  696  714  732  750  759  831  831  831  831'
  [../]

  # Table 2 in
  # Jun Cao, Michael A. Gharghouri, Philip Nash
  # Finite-element analysis and experimental validation of thermal residual stress and distortion in electron beam additive manufactured Ti-6Al-4V build plates
  # Journal of Materials Processing Technology, Volume 237, November 2016, Pages 409-419
  # https://www.sciencedirect.com/science/article/pii/S0924013616302126
  [./alpha_function]
    type = PiecewiseLinear
    x = '293    400     500    600    700    800    900   1000   1100    1200   1300  1400  1610   1950'
    y = '9.00e-06 9.16e-06 9.31e-06 9.46e-06 9.61e-06 9.76e-06 9.90e-06 1.01e-05 1.02e-05 1.04e-05 1.05e-05 1.06e-05 0.0 0.0'
  [../]

  # Yield stress as a function of temperature
  # Table 2 in
  # Jun Cao, Michael A. Gharghouri, Philip Nash
  # Finite-element analysis and experimental validation of thermal residual stress and distortion in electron beam additive manufactured Ti-6Al-4V build plates
  # Journal of Materials Processing Technology, Volume 237, November 2016, Pages 409-419
  # https://www.sciencedirect.com/science/article/pii/S0924013616302126
  [./yield_stress_function]
    type = PiecewiseLinear
    x = '293    400     500    600    700    800    900   1000   1100    1200   1300  1400    1950 19500'
    y = '8.50e8 7.20e8 6.80e8 6.30e8 5.90e8 5.40e8 4.90e8 4.50e8 4.00e8 3.60e8 3.15e8 2.68e8 2.00e7 2.00e7'
  [../]
  
  # x = plastic strain
  # y = increase in yield strength
  [./hardening_function]
    type = PiecewiseLinear
    x = '0.0 0.001 0.002 0.005 0.01'
    y = '0.0 90.9e6 118.1e6 143.3e6 170.0e6'    
  [../]
  
  [./timestep_fn]
    type = PiecewiseConstant
    x = '0.000000 0.000200 0.010000 0.010200 0.020000 0.020200 0.030000 0.030200 0.040000 0.040200 0.050000 0.050200 0.060000 0.060200 0.070000 0.070200 0.080000 0.080200 0.090000 0.090200 0.100000 0.100200 0.110000 0.110200 0.120000 0.120200 0.130000 0.130200 0.140000 0.140200 0.150000 0.150200 0.160000 0.160200 0.170000 0.170200 0.180000 0.180200 0.190000 0.190200'
    y = '1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5   1.2e-7   4.0e-5'
  [../]
  
  # element activation function
  [./level_set]
    type = ParsedFunction
    expression = 'if(z < 0.00014,1,0) * if(t > 0.0000003,1,0) + if(t > 0.099999,1,0)'
  [../]
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu' 
  petsc_options = '-snes_converged_reason'

  line_search = 'none'
  
  automatic_scaling = true

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25
    iteration_window = 15
    growth_factor = 1.02
    reject_large_step = True
    cutback_factor = 0.5
    timestep_limiting_postprocessor = timestep_pp
    dt = 1e-7
  []

  nl_max_its = 20
  l_max_its = 50
  
  start_time = 0.0
  end_time = 0.25
  dtmin = 1.0e-20
  
[]

[Postprocessors]
  [./temp]
    type = ElementAverageValue
    variable = temp
    block = 0
  [../]
  
  [./max_temperature]
    type = NodalExtremeValue
    variable = temp
    block = 0
  [../]

  [./stress_xx]
    type = ElementAverageValue
    variable = stress_xx
    block = 0
  [../]
  
  [./max_mechanical_strain_xx]
    type = ElementExtremeMaterialProperty
    value_type = max
    mat_prop = mechanical_strain_xx
    block = 0
  [../]
  
  [./max_mechanical_strain_yy]
    type = ElementExtremeMaterialProperty
    value_type = max
    mat_prop = mechanical_strain_yy
    block = 0
  [../]
  
  [./max_mechanical_strain_zz]
    type = ElementExtremeMaterialProperty
    value_type = max
    mat_prop = mechanical_strain_zz
    block = 0
  [../]
  
  [./timestep_pp]
    type = FunctionValuePostprocessor
    function = timestep_fn
  [../]

[]

# element activation user object
# element are moved to the subdomain 0
# where stress equilibrium and heat equation are solved
# at position and time where level_set becomes greater than 0.5 
[MeshModifiers]
  [activated_elem_uo] 
    type = ActivateElementsCoupled
    execute_on = timestep_begin
    activate_value = 0.5
    coupled_var = level_set
    active_subdomain_id = 0
    inactive_subdomain_id = 1
    expand_boundary_name = 'substrate_front'
  []
[]

[Outputs]

  [./out]
    type = Exodus
    interval = 1000
    #sync_only = true
    #sync_times = '0.000002 0.00001 0.00002 0.0001 0.001 0.0024 1.0'
  [../]

  csv = true
[]

