# A film of water under air freezes from a cold wall, with the ice less dense than the water
# (Lyu et al., J. Comput. Phys. 432 (2021) 110160, sec. 3). The expanding ice pushes the water,
# and the air above it, upwards. While the water is deep against the thermal boundary layer, the
# ice front lies at
#   h_s = 2 delta sqrt(a_s t),
# with delta the root of Lyu's eq. (30),
#   exp(-delta^2) / erf(delta)
#     - phi sqrt(a) exp(-[sqrt(a) delta (1 - r)]^2) / erfc[sqrt(a) delta (1 - r)]
#     = delta sqrt(pi) / Ste,
# where a = a_s / a_l is the ratio of the thermal diffusivities, phi = k_l (T_i - T_m) / (k_s (T_m -
# T_w)), r = 1 - rho_s / rho_l and Ste = cp_s (T_m - T_w) / L. The water above the front moves at
# r dh_s/dt. Once all the water has frozen, the ice is h_0 rho_l / rho_s high. The problem is
# one-dimensional, computed on a column two cells wide.

# Water, ice and air (Lyu et al., Table 1); the conductivities are rho cp times the diffusivity
rho_l = 1000
cp_l = 4210
k_l = ${fparse rho_l * cp_l * 1.32e-7}
mu_l = 1.7e-3
rho_s = 917
cp_s = 2030
k_s = ${fparse rho_s * cp_s * 1.18e-6}
rho_g = 1.29
cp_g = 1000
k_g = ${fparse rho_g * cp_g * 2e-5}
mu_g = 1.7e-5
L = 334e3
T_m = 273.15
T_w = ${fparse T_m - 20}
T_i = ${fparse T_m + 20}

h_0 = 1e-3
a_s = ${fparse k_s / (rho_s * cp_s)}
# Root of the equation above for these properties, to 12 digits
delta = 0.203004067523

# Column of 200 cells over twice the water depth, two cells wide
dy = ${fparse 2 * h_0 / 200}
width = ${fparse 2 * dy}

# Brinkman penalty holding the ice at rest: its coefficient rho alpha (1 - f_l) / (eta dt) is
# 1 / eta times the time term of the momentum equation in a solid cell
eta = 1e-3

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${width}'
    dy = '${fparse 2 * h_0}'
    ix = '2'
    iy = '200'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system alpha_system energy_system'
  previous_nl_solution_required = true
[]

[Physics]
  [NavierStokes]
    [ConservativeDiffuseInterfaceFlowSegregated]
      [flow]
        velocity_variable = 'vel_x vel_y'
        pressure_variable = 'pressure'

        compressibility = 'incompressible'
        density = 'rho_mixture'
        dynamic_viscosity = 'mu_mixture'
        gravity = '0 0 0'
        volume_fraction_functor = 'alpha'

        initial_velocity = '0 0 0'
        initial_pressure = 0

        wall_boundaries = 'left right bottom'
        momentum_wall_types = 'symmetry symmetry noslip'

        outlet_boundaries = 'top'
        momentum_outlet_types = 'fixed-pressure'
        pressure_functors = '0'

        orthogonality_correction = false
        momentum_two_term_bc_expansion = false
        pressure_two_term_bc_expansion = false
        momentum_advection_interpolation = 'upwind'
      []
    []
    [ConservativeDiffuseInterfaceVOFSegregated]
      [vof]
        coupled_flow_physics = 'flow'
        passive_scalar_names = 'alpha'
        initial_scalar_variables = 'alpha_init'
        system_names = 'alpha_system'

        liquid_density_name = 'rho_l'
        gas_density_name = 'rho_g'
        liquid_dynamic_viscosity_name = 'mu_l'
        gas_dynamic_viscosity_name = 'mu_g'

        passive_scalar_advection_interpolation = 'upwind'
        compression_factor = 0.1
        interface_normal_functor = 'flow_interface_unit_normal_face'

        n_alpha_corrections = 1
        n_limiter_iterations = 6
      []
    []
    [ConservativeDiffuseInterfaceHeatTransferSegregated]
      [heat]
        coupled_vof_physics = 'vof'
        conserved_enthalpy_variable = 'rho_h'
        temperature_variable = 'T'
        liquid_fluid_properties = 'fp_water'
        gas_fluid_properties = 'fp_air'
        solid_fluid_properties = 'fp_ice'
        solid_density = 'rho_s'
        system_names = 'energy_system'
        thermodynamic_pressure = 0

        initial_temperature = '${T_i}'
        latent_heat = ${L}
        solidus_temperature = ${T_m}
        liquidus_temperature = ${T_m}

        energy_wall_types = 'heatflux heatflux fixed-temperature'
        energy_wall_functors = '0 0 ${T_w}'
      []
    []
  []
[]

[FluidProperties]
  [fp_water]
    type = TemperaturePressureFunctionFluidProperties
    rho = ${rho_l}
    mu = ${mu_l}
    k = ${k_l}
    cv = ${cp_l}
    T_ref = ${T_m}
  []
  [fp_ice]
    type = TemperaturePressureFunctionFluidProperties
    rho = ${rho_s}
    mu = ${mu_l}
    k = ${k_s}
    cv = ${cp_s}
    T_ref = ${T_m}
  []
  [fp_air]
    type = TemperaturePressureFunctionFluidProperties
    rho = ${rho_g}
    mu = ${mu_g}
    k = ${k_g}
    cv = ${cp_g}
    T_ref = ${T_m}
  []
[]

[Functions]
  [alpha_init]
    type = ParsedFunction
    expression = 'if(y < ${h_0}, 1, 0)'
  []
  [h_s_exact]
    type = ParsedFunction
    expression = '2 * ${delta} * sqrt(${a_s} * t)'
  []
[]

[Postprocessors]
  [dt_pp]
    type = TimestepSize
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[FunctorMaterials]
  [constants]
    type = GenericFunctorMaterial
    prop_names = 'rho_l rho_s rho_g mu_l mu_g'
    prop_values = '${rho_l} ${rho_s} ${rho_g} ${mu_l} ${mu_g}'
  []
  [solid_fraction]
    type = ParsedFunctorMaterial
    property_name = 'solid_fraction'
    functor_names = 'alpha liquid_fraction'
    expression = 'alpha * (1 - liquid_fraction)'
  []
  [penalty]
    type = ParsedFunctorMaterial
    property_name = 'penalty_coefficient'
    functor_names = 'rho_mixture solid_fraction dt_pp'
    expression = 'rho_mixture * solid_fraction / (${eta} * dt_pp)'
  []
  # Volume source of the freezing water: a cell keeps the mass of the ice it forms and expels the
  # rest as water, so div(u) = alpha (1 - rho_s / rho_l) df_s/dt (Lyu et al., eq. 11)
  [liquid_fraction_change]
    type = FunctorChangeFunctorMaterial
    functor = 'liquid_fraction'
    change_over = 'time_step'
    take_absolute_value = false
    prop_name = 'liquid_fraction_change'
  []
  [expansion_rate]
    type = ParsedFunctorMaterial
    property_name = 'expansion_rate'
    functor_names = 'liquid_fraction_change dt_pp'
    functor_symbols = 'df_l dt'
    expression = '-(1 - ${rho_s} / ${rho_l}) * df_l / dt'
  []
  [volume_source]
    type = ParsedFunctorMaterial
    property_name = 'volume_source'
    functor_names = 'alpha expansion_rate'
    expression = 'alpha * expansion_rate'
  []
[]

[LinearFVKernels]
  [u_penalty]
    type = LinearFVReaction
    variable = vel_x
    coeff = 'penalty_coefficient'
  []
  [v_penalty]
    type = LinearFVReaction
    variable = vel_y
    coeff = 'penalty_coefficient'
  []
  [pressure_volume_source]
    type = LinearFVSource
    variable = pressure
    source_density = 'volume_source'
  []
  # The volume the freezing water gains is water: d(alpha)/dt + div(alpha u) = alpha div(u), taken
  # as the divergence of the same face flux that advects alpha so that alpha stays exactly 1 in
  # the water whatever the state of the outer iteration
  [alpha_volume_source]
    type = LinearFVDivergence
    variable = alpha
    face_flux = 'lagged_corrected_face_phi'
    face_flux_is_integrated = false
    force_boundary_execution = true
  []
[]

[Executioner]
  type = ReducedPressurePIMPLE
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'

  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  active_scalar_systems = 'alpha_system'
  energy_system = 'energy_system'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  active_scalar_equation_relaxation = '1.0'
  energy_equation_relaxation = 1.0

  num_iterations = 20
  num_piso_iterations = 0
  continue_on_max_its = true
  print_fields = false

  momentum_absolute_tolerance = 1e-10
  pressure_absolute_tolerance = 1e-10
  active_scalar_absolute_tolerance = '1e-10'
  energy_absolute_tolerance = 1e-6

  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  active_scalar_l_abs_tol = 1e-12
  energy_l_abs_tol = 1e-12

  momentum_l_tol = 0
  pressure_l_tol = 0
  active_scalar_l_tol = 0
  energy_l_tol = 0

  momentum_petsc_options_iname = '-pc_type'
  momentum_petsc_options_value = 'lu'
  pressure_petsc_options_iname = '-pc_type'
  pressure_petsc_options_value = 'lu'
  active_scalar_petsc_options_iname = '-pc_type'
  active_scalar_petsc_options_value = 'lu'
  energy_petsc_options_iname = '-pc_type'
  energy_petsc_options_value = 'lu'

  dt = 2e-3
  # The front has crossed 25 cells; the water above is still deep against its thermal boundary
  # layer, 2 sqrt(a_l t) = 0.4 mm
  end_time = 0.3
[]

[Postprocessors]
  # Ice and water depths per unit wall area
  [h_s]
    type = ElementIntegralFunctorPostprocessor
    functor = 'solid_fraction'
    prefactor = ${fparse 1 / width}
  []
  [h_s_exact]
    type = FunctionValuePostprocessor
    function = 'h_s_exact'
  []
  [h_s_error]
    type = ParsedPostprocessor
    expression = 'h_s / h_s_exact - 1'
    pp_names = 'h_s h_s_exact'
  []
  # Water and ice together: the ice formed so far has expelled r h_s of water
  [h_liquid_phase]
    type = ElementIntegralFunctorPostprocessor
    functor = 'alpha'
    prefactor = ${fparse 1 / width}
  []
  [h_liquid_phase_error]
    type = ParsedPostprocessor
    expression = 'h_liquid_phase / (${h_0} + (1 - ${rho_s} / ${rho_l}) * h_s) - 1'
    pp_names = 'h_liquid_phase h_s'
  []
  [alpha_max]
    type = ElementExtremeValue
    variable = 'alpha'
    value_type = max
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
