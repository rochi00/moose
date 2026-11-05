# A liquid drop on a wall relaxes from a hemisphere to the spherical cap of its contact angle,
# in axisymmetric coordinates without gravity. For a cap of volume V and contact angle theta, the
# radius, height and base radius are
#   R = (3 V / (pi (2 + cos theta) (1 - cos theta)^2))^(1/3), H = R (1 - cos theta),
#   a = R sin theta.
# The height is sampled as the liquid volume of the first column of cells over its volume, the base
# radius from the wetted area of the wall. rho_l = rho_g, and the viscosity damps the capillary
# oscillation in a few periods: the Ohnesorge number mu / sqrt(rho sigma R_0) is 0.14.

theta = 60
R_0 = 0.5
sigma = 1
mu = 1
h = ${fparse 1 / 32}
interface_width = ${fparse 4 * h}

cos_theta = ${fparse cos(theta * pi / 180)}
sin_theta = ${fparse sin(theta * pi / 180)}
# Cap radius per unit cube root of the liquid volume, which the diffuse hemisphere carries more of
# than a sharp one
R_cap_per_V3 = ${fparse (3 / (pi * (2 + cos_theta) * (1 - cos_theta)^2))^(1 / 3)}

[Mesh]
  coord_type = RZ
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1'
    dy = '1'
    ix = '32'
    iy = '32'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system alpha_system'
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
        surface_tension_coefficient = 'sigma'
        contact_angle = ${theta}
        contact_angle_boundaries = 'bottom'

        initial_velocity = '0 0 0'
        initial_pressure = '0'

        wall_boundaries = 'left right bottom'
        momentum_wall_types = 'symmetry noslip symmetry'

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

        liquid_density_name = 'rho'
        gas_density_name = 'rho'
        liquid_dynamic_viscosity_name = 'mu'
        gas_dynamic_viscosity_name = 'mu'

        passive_scalar_advection_interpolation = 'upwind'
        compression_factor = 0.1
        interface_normal_functor = 'flow_interface_unit_normal_face'

        n_alpha_corrections = 1
        n_limiter_iterations = 6
      []
    []
  []
[]

[FunctorMaterials]
  [constants]
    type = GenericFunctorMaterial
    prop_names = 'rho mu sigma'
    prop_values = '1 ${mu} ${sigma}'
  []
[]

[Functions]
  [alpha_init]
    type = ParsedFunction
    expression = '0.5 * (1 + tanh((${R_0} - sqrt(x^2 + y^2)) / (sqrt(2) * ${interface_width})))'
  []
[]

[AuxVariables]
  [axis_liquid]
    type = MooseVariableFVReal
  []
  [axis_column]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  # The first column of cells next to the axis
  [axis_liquid]
    type = ParsedAux
    variable = axis_liquid
    coupled_variables = 'alpha'
    expression = 'alpha * if(x < ${h}, 1, 0)'
    use_xyzt = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [axis_column]
    type = ParsedAux
    variable = axis_column
    expression = 'if(x < ${h}, 1, 0)'
    use_xyzt = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = ReducedPressurePIMPLE
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'

  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  active_scalar_systems = 'alpha_system'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 1.0
  active_scalar_equation_relaxation = '1.0'

  num_iterations = 4
  num_piso_iterations = 0
  continue_on_max_its = true
  print_fields = false

  momentum_absolute_tolerance = 1e-6
  pressure_absolute_tolerance = 1e-6
  active_scalar_absolute_tolerance = '1e-6'

  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  active_scalar_l_abs_tol = 1e-12

  momentum_l_tol = 0
  pressure_l_tol = 0
  active_scalar_l_tol = 0

  momentum_petsc_options_iname = '-pc_type'
  momentum_petsc_options_value = 'lu'
  pressure_petsc_options_iname = '-pc_type'
  pressure_petsc_options_value = 'lu'
  active_scalar_petsc_options_iname = '-pc_type'
  active_scalar_petsc_options_value = 'lu'

  startup_flux_corrections = 2
  volume_fraction_subcycles = 1

  # Below the capillary time step sqrt(rho h^3 / (2 pi sigma)) = 2.2e-3 of the explicit surface force
  dt = 2e-3
  end_time = 3
[]

[Postprocessors]
  [liquid_volume]
    type = ElementIntegralVariablePostprocessor
    variable = 'alpha'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [axis_liquid_volume]
    type = ElementIntegralVariablePostprocessor
    variable = 'axis_liquid'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [axis_column_volume]
    type = ElementIntegralVariablePostprocessor
    variable = 'axis_column'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [height]
    type = ParsedPostprocessor
    pp_names = 'axis_liquid_volume axis_column_volume'
    expression = 'axis_liquid_volume / axis_column_volume'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [height_error]
    type = ParsedPostprocessor
    pp_names = 'height liquid_volume'
    expression = 'height / (${R_cap_per_V3} * liquid_volume^(1/3) * (1 - ${cos_theta})) - 1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [wetted_area]
    type = SideIntegralVariablePostprocessor
    variable = 'alpha'
    boundary = 'bottom'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [base_radius]
    type = ParsedPostprocessor
    pp_names = 'wetted_area'
    expression = 'sqrt(wetted_area / ${fparse pi})'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [base_radius_error]
    type = ParsedPostprocessor
    pp_names = 'base_radius liquid_volume'
    expression = 'base_radius / (${R_cap_per_V3} * liquid_volume^(1/3) * ${sin_theta}) - 1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
