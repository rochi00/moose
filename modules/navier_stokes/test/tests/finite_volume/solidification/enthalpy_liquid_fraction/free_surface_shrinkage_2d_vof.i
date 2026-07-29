rho_gas = 0.01
rho_solid = 1.2
rho_liquid = 1
k_gas = 0.01
k_solid = 1
k_liquid = 0.05
cp_gas = 1
cp_solid = 1
cp_liquid = 1
L = 0.53
T_solidus = 1
T_liquidus = 1.01
T_gas = 1.0
T_cold = 0.8
T_solid = 0.9
T_liquid = 1.2
solid_liquid_interface_y = 0.15
material_interface_y = 0.70
material_interface_amplitude = 0.03
domain_width = 1
mu_gas = 1
mu = 1
DarcyConstantlarge = 1e6
DarcyConstantsmall = 1e-3
Nx = 32
Ny = 80
c_alpha = 0.0

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${domain_width}
    ymin = 0
    ymax = 1
    nx = ${Nx}
    ny = ${Ny}
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system alpha_system energy_system'
  previous_nl_solution_required = true
[]

[Functions]
  [initial_alpha]
    type = ParsedFunction
    expression = 'if(y < material_interface_y + material_interface_amplitude * cos(2 * pi * x / domain_width), 1, 0)'
    symbol_names = 'material_interface_y material_interface_amplitude domain_width'
    symbol_values = '${material_interface_y} ${material_interface_amplitude} ${domain_width}'
  []
  [initial_fl]
    type = ParsedFunction
    expression = 'if(y < material_interface_y + material_interface_amplitude * cos(2 * pi * x / domain_width), if(y < solid_liquid_interface_y, 0, 1), 0)'
    symbol_names = 'material_interface_y material_interface_amplitude domain_width solid_liquid_interface_y'
    symbol_values = '${material_interface_y} ${material_interface_amplitude} ${domain_width} ${solid_liquid_interface_y}'
  []
  [initial_T]
    type = ParsedFunction
    expression = 'if(y < material_interface_y + material_interface_amplitude * cos(2 * pi * x / domain_width), if(y < solid_liquid_interface_y, T_solid, T_liquid), T_gas)'
    symbol_names = 'material_interface_y material_interface_amplitude domain_width solid_liquid_interface_y T_solid T_liquid T_gas'
    symbol_values = '${material_interface_y} ${material_interface_amplitude} ${domain_width} ${solid_liquid_interface_y} ${T_solid} ${T_liquid} ${T_gas}'
  []
  [initial_rho_h]
    type = ParsedFunction
    expression = 'if(y < material_interface_y + material_interface_amplitude * cos(2 * pi * x / domain_width), if(y < solid_liquid_interface_y, rho_solid * cp_solid * T_solid, rho_liquid * cp_liquid * T_liquid), rho_gas * cp_gas * T_gas)'
    symbol_names = 'material_interface_y material_interface_amplitude domain_width solid_liquid_interface_y rho_solid rho_liquid rho_gas cp_solid cp_liquid cp_gas T_solid T_liquid T_gas'
    symbol_values = '${material_interface_y} ${material_interface_amplitude} ${domain_width} ${solid_liquid_interface_y} ${rho_solid} ${rho_liquid} ${rho_gas} ${cp_solid} ${cp_liquid} ${cp_gas} ${T_solid} ${T_liquid} ${T_gas}'
  []
  [bottom_rho_h]
    type = ParsedFunction
    expression = 'rho_solid * cp_solid * T_cold'
    symbol_names = 'rho_solid cp_solid T_cold'
    symbol_values = '${rho_solid} ${cp_solid} ${T_cold}'
  []
  [pressure_init]
    type = ConstantFunction
    value = 0
  []
[]

[AuxVariables]
  [fl]
    type = MooseVariableFVReal
  []
  [density]
    type = MooseVariableFVReal
  []
  [phase_change_divergence_out]
    type = MooseVariableFVReal
  []
  [alpha_liquid]
    type = MooseVariableFVReal
  []
  [alpha_solid]
    type = MooseVariableFVReal
  []
[]

[ICs]
  [fl_ic]
    type = FunctionIC
    variable = fl
    function = initial_fl
  []
[]

[UserObjects]
  [liquid_fraction_corrector]
    type = EnthalpyLiquidFractionCorrector
    liquid_fraction_variable = fl
    temperature = temperature
    cp = cp_pcm
    L = latent_heat_material
    material_fraction = alpha
    min_material_fraction = 0.5
    T_liquidus = ${T_liquidus}
    T_solidus = ${T_solidus}
    relaxation = 0.5
    max_liquid_fraction_change = 0.005
    phase_change_rate_relaxation = 0.1
    execute_on = 'INITIAL'
  []
[]

[Physics]
  [NavierStokes]
    [ConservativeSharpInterfaceFlowSegregated]
      [flow]
        velocity_variable = 'vel_x vel_y'
        pressure_variable = 'pressure'

        compressibility = 'incompressible'
        density = 'rho_mixture'
        dynamic_viscosity = 'mu_mixture'
        gravity = '0 0 0'
        volume_fraction_functor = 'alpha'

        initial_velocity = '0 0 0'
        initial_pressure = 'pressure_init'

        wall_boundaries = 'left right bottom'
        momentum_wall_types = 'noslip noslip noslip'

        outlet_boundaries = 'top'
        momentum_outlet_types = 'fixed-pressure'
        pressure_functors = '0'

        orthogonality_correction = false
        momentum_two_term_bc_expansion = false
        pressure_two_term_bc_expansion = false
        momentum_advection_interpolation = 'upwind'
      []
    []
    [ConservativeSharpInterfaceVOFSegregated]
      [vof]
        coupled_flow_physics = 'flow'
        passive_scalar_names = 'alpha'
        initial_scalar_variables = 'initial_alpha'
        system_names = 'alpha_system'
        volume_fraction_outlet_type = 'inlet-outlet'

        liquid_density_name = 'rho_pcm'
        gas_density_name = 'rho_gas_mat'
        liquid_specific_heat_name = 'cp_pcm'
        gas_specific_heat_name = 'cp_gas_mat'
        rho_cp_phi_name = 'rho_cp_phi'
        liquid_dynamic_viscosity_name = 'mu_pcm'
        gas_dynamic_viscosity_name = 'mu_gas_mat'

        passive_scalar_advection_interpolation = 'upwind'
        compression_factor = '${c_alpha}'
        interface_normal_functor = 'flow_interface_unit_normal_face'
        phase_change_divergence = 'phase_change_divergence'
        conserved_enthalpy_variable = 'rho_h'
        conserved_enthalpy_temperature = 'temperature'
        conserved_enthalpy_backflow_temperature = 'T_gas_mat'

        n_alpha_corrections = 1
        n_limiter_iterations = 6
      []
    []
    [FluidHeatTransferSegregated]
      [heat]
        coupled_flow_physics = 'flow'
        system_names = 'energy_system'
        fluid_temperature_variable = 'temperature'
        solve_for_conserved_enthalpy = true
        fluid_conserved_enthalpy_variable = 'rho_h'

        thermal_conductivity = 'k_mixture'
        specific_heat = 'cp_mixture'
        energy_mass_heat_capacity_face_flux = 'rho_cp_phi'
        energy_mass_heat_capacity_face_flux_is_integrated = true
        use_vof_consistent_energy_advection = true

        initial_conserved_enthalpy = 'initial_rho_h'

        energy_wall_boundaries = 'bottom'
        energy_wall_types = 'fixed-temperature'
        energy_wall_functors = 'bottom_rho_h'

        use_nonorthogonal_correction = false
        energy_two_term_bc_expansion = false
        energy_advection_interpolation = 'upwind'
      []
    []
  []
[]

[AuxKernels]
  [rho_out]
    type = FunctorAux
    variable = density
    functor = rho_mixture
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pc_div_out]
    type = FunctorAux
    variable = phase_change_divergence_out
    functor = phase_change_divergence
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [alpha_liquid_out]
    type = FunctorAux
    variable = alpha_liquid
    functor = alpha_liquid_material
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [alpha_solid_out]
    type = FunctorAux
    variable = alpha_solid
    functor = alpha_solid_material
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[LinearFVKernels]
  [u_friction_darcy]
    type = LinearFVReaction
    variable = vel_x
    coeff = darcy_coef_friction
  []
  [v_friction_darcy]
    type = LinearFVReaction
    variable = vel_y
    coeff = darcy_coef_friction
  []
  [phase_change_divergence]
    type = LinearFVSource
    variable = pressure
    source_density = phase_change_divergence
    scaling_factor = 1
  []
  [latent_time]
    type = LinearFVEnthalpyLiquidFractionLatentHeatTimeDerivative
    variable = rho_h
    L = latent_heat_material
    density = rho_mixture
    liquid_fraction = fl
  []
  [latent_advection]
    type = LinearFVEnthalpyLiquidFractionLatentHeatAdvection
    variable = rho_h
    L = latent_heat_material
    liquid_fraction = fl
    mass_flux_functor = rho_phi
    mass_flux_is_integrated = true
  []
[]

[FunctorMaterials]
  [constants]
    type = GenericFunctorMaterial
    prop_names = 'rho_gas_mat rho_solid_mat rho_liquid_mat cp_gas_mat cp_solid_mat cp_liquid_mat k_gas_mat k_solid_mat k_liquid_mat mu_gas_mat mu_liquid_mat T_gas_mat'
    prop_values = '${rho_gas} ${rho_solid} ${rho_liquid} ${cp_gas} ${cp_solid} ${cp_liquid} ${k_gas} ${k_solid} ${k_liquid} ${mu_gas} ${mu} ${T_gas}'
  []
  [rho_pcm_mat]
    type = ParsedFunctorMaterial
    property_name = rho_pcm
    functor_names = 'fl rho_solid_mat rho_liquid_mat'
    functor_symbols = 'fl rho_solid rho_liquid'
    expression = 'rho_solid * (1 - min(max(fl, 0), 1)) + rho_liquid * min(max(fl, 0), 1)'
  []
  [cp_pcm_mat]
    type = ParsedFunctorMaterial
    property_name = cp_pcm
    functor_names = 'fl cp_solid_mat cp_liquid_mat'
    functor_symbols = 'fl cp_solid cp_liquid'
    expression = 'cp_solid * (1 - min(max(fl, 0), 1)) + cp_liquid * min(max(fl, 0), 1)'
  []
  [mu_pcm_mat]
    type = ParsedFunctorMaterial
    property_name = mu_pcm
    functor_names = 'fl mu_liquid_mat'
    functor_symbols = 'fl mu_liquid'
    expression = 'mu_liquid'
  []
  [rho_cp_props]
    type = ParsedFunctorMaterial
    property_name = rho_cp
    functor_names = 'alpha fl rho_gas_mat cp_gas_mat rho_solid_mat cp_solid_mat rho_liquid_mat cp_liquid_mat'
    functor_symbols = 'alpha fl rho_gas cp_gas rho_solid cp_solid rho_liquid cp_liquid'
    expression = 'min(max(alpha, 0), 1) * (rho_solid * cp_solid * (1 - min(max(fl, 0), 1)) + rho_liquid * cp_liquid * min(max(fl, 0), 1)) + (1 - min(max(alpha, 0), 1)) * rho_gas * cp_gas'
  []
  [cp_mixture_mat]
    type = ParsedFunctorMaterial
    property_name = cp_mixture
    functor_names = 'rho_cp rho_mixture'
    functor_symbols = 'rho_cp rho'
    expression = 'rho_cp / rho'
  []
  [k_mixture_mat]
    type = ParsedFunctorMaterial
    property_name = k_mixture
    functor_names = 'alpha fl k_gas_mat k_solid_mat k_liquid_mat'
    functor_symbols = 'alpha fl k_gas k_solid k_liquid'
    expression = 'min(max(alpha, 0), 1) * (k_solid * (1 - min(max(fl, 0), 1)) + k_liquid * min(max(fl, 0), 1)) + (1 - min(max(alpha, 0), 1)) * k_gas'
  []
  [latent_heat]
    type = ParsedFunctorMaterial
    property_name = latent_heat_material
    functor_names = 'alpha'
    functor_symbols = 'alpha'
    expression = 'min(max(alpha, 0), 1) * ${L}'
  []
  [phase_change_volume_change]
    type = EnthalpyLiquidFractionPhaseChangeDivergenceFunctorMaterial
    liquid_fraction = fl
    density = rho_mixture
    rho_liquid = ${rho_liquid}
    rho_solid = ${rho_solid}
    material_fraction = alpha
    liquid_fraction_corrector = liquid_fraction_corrector
  []
  [darcy_coeff_friction]
    type = ParsedFunctorMaterial
    property_name = darcy_coef_friction
    functor_names = 'alpha fl'
    functor_symbols = 'alpha fl'
    expression = 'min(max(alpha, 0), 1) * ${DarcyConstantlarge} * (1 - min(max(fl, 0), 1))^2 / (min(max(fl, 0), 1)^3 + ${DarcyConstantsmall})'
  []
  [alpha_liquid_fraction]
    type = ParsedFunctorMaterial
    property_name = alpha_liquid_material
    functor_names = 'alpha fl'
    functor_symbols = 'alpha fl'
    expression = 'min(max(alpha, 0), 1) * min(max(fl, 0), 1)'
  []
  [alpha_solid_fraction]
    type = ParsedFunctorMaterial
    property_name = alpha_solid_material
    functor_names = 'alpha fl'
    functor_symbols = 'alpha fl'
    expression = 'min(max(alpha, 0), 1) * (1 - min(max(fl, 0), 1))'
  []
[]

[Executioner]
  type = ReducedPressurePIMPLE
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'

  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  energy_system = 'energy_system'
  active_scalar_systems = 'alpha_system'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.05
  energy_equation_relaxation = 0.9
  active_scalar_equation_relaxation = '1'

  liquid_fraction_corrector = liquid_fraction_corrector
  min_temperature_correctors = 1
  max_temperature_correctors = 20
  liquid_fraction_tolerance = 1e-4

  num_iterations = 30
  num_piso_iterations = 0
  continue_on_max_its = true
  print_fields = false

  momentum_absolute_tolerance = 1e-10
  pressure_absolute_tolerance = 1e-10
  energy_absolute_tolerance = 1e-10
  active_scalar_absolute_tolerance = '1e-10'

  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  energy_l_abs_tol = 1e-12
  active_scalar_l_abs_tol = 1e-12

  momentum_l_tol = 0
  pressure_l_tol = 0
  energy_l_tol = 0
  active_scalar_l_tol = 0

  momentum_petsc_options_iname = '-pc_type'
  momentum_petsc_options_value = 'lu'
  pressure_petsc_options_iname = '-pc_type'
  pressure_petsc_options_value = 'lu'
  energy_petsc_options_iname = '-pc_type'
  energy_petsc_options_value = 'lu'
  active_scalar_petsc_options_iname = '-pc_type'
  active_scalar_petsc_options_value = 'lu'

  startup_pressure_initialization = 'projection-only'
  startup_flux_corrections = 2
  volume_fraction_subcycles = 2

  dt = 1e-4
  num_steps = 10
[]

[Postprocessors]
  [alpha_min]
    type = ElementExtremeValue
    variable = alpha
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [alpha_max]
    type = ElementExtremeValue
    variable = alpha
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_alpha]
    type = ElementIntegralVariablePostprocessor
    variable = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_liquid]
    type = ElementIntegralVariablePostprocessor
    variable = alpha_liquid
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_solid]
    type = ElementIntegralVariablePostprocessor
    variable = alpha_solid
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ave_T]
    type = ElementAverageValue
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ave_pc_div]
    type = ElementAverageValue
    variable = phase_change_divergence_out
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ave_density]
    type = ElementAverageValue
    variable = density
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  execute_on = 'INITIAL TIMESTEP_END'
  csv = true
  exodus = true
[]
