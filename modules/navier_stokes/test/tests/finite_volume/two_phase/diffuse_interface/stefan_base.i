# Shared setup of the Stefan melting tests, included after the including input has defined
# T_s and T_l and the functions T_exact and melt_exact.
#
# A semi-infinite solid at T_i, at rest and filling the domain (alpha = 1), melts from a wall
# raised to T_w at x = 0 at t = 0. rho = cp = k = L = 1, so the Stefan numbers are 1 on the liquid
# side and 0.5 on the solid side. The far wall at x = 1 holds the exact (solid) temperature. The
# similarity solution holds from t = 0, where T_exact is T_i in the interior, so it also sets the
# initial condition; the run is compared with it at t_1. h = T - 300 exactly: cv is constant, the
# reference temperature is 300 and the thermodynamic pressure zero.

T_w = 301
T_i = 299.5
L = 1
t_1 = 0.032

[Mesh]
  # The problem is one-dimensional. Two cells across, as with one the cell gradient reconstruction
  # couples the top and bottom boundaries of every cell and is singular.
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1'
    dy = '1'
    ix = '100'
    iy = '2'
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
    [ConservativeDiffuseInterfaceVOFSegregated]
      [vof]
        coupled_flow_physics = 'flow'
        passive_scalar_names = 'alpha'
        initial_scalar_variables = '1'
        system_names = 'alpha_system'

        liquid_density_name = 'rho'
        gas_density_name = 'rho'
        liquid_dynamic_viscosity_name = 'mu'
        gas_dynamic_viscosity_name = 'mu'

        passive_scalar_advection_interpolation = 'upwind'
        compression_factor = 0
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
        liquid_fluid_properties = 'fp'
        gas_fluid_properties = 'fp'
        system_names = 'energy_system'
        thermodynamic_pressure = 0

        initial_temperature = 'T_exact'
        latent_heat = ${L}
        solidus_temperature = ${T_s}
        liquidus_temperature = ${T_l}

        energy_wall_types = 'fixed-temperature fixed-temperature heatflux'
        energy_wall_functors = '${T_w} T_exact 0'
      []
    []
  []
[]

[FluidProperties]
  [fp]
    type = TemperaturePressureFunctionFluidProperties
    rho = 1
    mu = 1
    k = 1
    cv = 1
    T_ref = 300
  []
[]

[FunctorMaterials]
  [constants]
    type = GenericFunctorMaterial
    prop_names = 'rho mu'
    prop_values = '1 1'
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

  # Each step is iterated until rho*h and the recovered temperature and liquid fraction agree
  num_iterations = 20
  num_piso_iterations = 0
  continue_on_max_its = true
  print_fields = false

  momentum_absolute_tolerance = 1e-10
  pressure_absolute_tolerance = 1e-10
  active_scalar_absolute_tolerance = '1e-10'
  energy_absolute_tolerance = 1e-10

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

  end_time = ${t_1}
  dt = 2e-4
[]

[Postprocessors]
  # L2 error of the cell temperatures against the exact profile; the domain is one unit high, so
  # integrals over it are per unit wall area
  [T_error]
    type = ElementL2FunctorError
    approximate = 'T'
    exact = 'T_exact'
  []
  # Melted volume per unit wall area, and its relative error
  [melt]
    type = ElementIntegralVariablePostprocessor
    variable = 'liquid_fraction'
  []
  [melt_exact]
    type = FunctionValuePostprocessor
    function = 'melt_exact'
  []
  [melt_error]
    type = ParsedPostprocessor
    expression = '(melt - melt_exact) / melt_exact'
    pp_names = 'melt melt_exact'
  []
[]

[Outputs]
  # Only the csv output is restricted to the end: a common execute_on would also apply to the
  # checkpoints of a recover run
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
