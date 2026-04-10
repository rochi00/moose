# Single-phase VOF static droplet: Laplace pressure validation
# One velocity field, mixture properties, algebraic MULES advection for alpha
# Phase 1 of CLSVOF implementation plan

rho1 = 10.0
rho2 = 1.0
mu_1 = 1e-3
mu_2 = 1e-3
surface_tension = 0.072

# Square domain and centered droplet
L = 1.0
R = ${fparse L / 4.0}
xc = ${fparse L / 2.0}
yc = ${fparse L / 2.0}
eps = ${fparse 4.0 * L / nxy}

# Keep boundary probes away from interface
inside_probe_half = ${fparse R / 4.0}
outside_probe_half = ${fparse L / 10.0}

# Mesh resolution
nxy = 100

advected_interp_method = 'upwind'
limiter_method = 'quick'

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${L}'
    dy = '${L}'
    ix = '${nxy}'
    iy = '${nxy}'
  []
  [inside_probe]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 1
    bottom_left = '${fparse xc - inside_probe_half} ${fparse yc - inside_probe_half} 0'
    top_right = '${fparse xc + inside_probe_half} ${fparse yc + inside_probe_half} 0'
  []
  [outside_probe]
    type = SubdomainBoundingBoxGenerator
    input = inside_probe
    block_id = 2
    bottom_left = '${fparse 0.0} ${fparse 0.0} 0'
    top_right = '${fparse outside_probe_half} ${fparse outside_probe_half} 0'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system alpha_1_system'
  previous_nl_solution_required = true
[]

[UserObjects]
  [rc]
    type = RhieChowMassFluxMultiPhase
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = rho_mix
    p_diffusion_kernel = p_diffusion
    check_executioner = false
    alpha = 1.0
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_condition = 0.0
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_condition = 0.0
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0.0
  []
  [alpha_1]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_1_system
  []
[]

[FunctorMaterials]
  [mixture]
    type = WCNSLinearFVMixtureFunctorMaterial
    phase_1_names = '${rho1} ${mu_1}'
    phase_2_names = '${rho2} ${mu_2}'
    prop_names = 'rho_mix mu_mix'
    phase_1_fraction = alpha_1
    limit_phase_fraction = true
  []
[]

[LinearFVKernels]
  # Momentum: d(rho_mix*u)/dt + div(rho_mix*u*u) = -grad(p) + div(mu_mix*grad(u)) + F_sigma
  [u_time]
    type = LinearFVTimeDerivative
    variable = vel_x
    factor = rho_mix
  []
  [v_time]
    type = LinearFVTimeDerivative
    variable = vel_y
    factor = rho_mix
  []
  [u_advection_stress]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    mu = mu_mix
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    rhie_chow_user_object = 'rc'
    alpha = 1.0
    use_nonorthogonal_correction = false
  []
  [v_advection_stress]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    mu = mu_mix
    u = vel_x
    v = vel_y
    momentum_component = 'y'
    rhie_chow_user_object = 'rc'
    alpha = 1.0
    use_nonorthogonal_correction = false
  []
  [u_pressure]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_x
    pressure = pressure
    alpha = 1.0
    momentum_component = 'x'
  []
  [v_pressure]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_y
    pressure = pressure
    alpha = 1.0
    momentum_component = 'y'
  []
  [u_surface_tension]
    type = LinearFVMomentumSurfaceTensionForce
    variable = vel_x
    sigma = ${surface_tension}
    alpha = alpha_1
    momentum_component = 'x'
  []
  [v_surface_tension]
    type = LinearFVMomentumSurfaceTensionForce
    variable = vel_y
    sigma = ${surface_tension}
    alpha = alpha_1
    momentum_component = 'y'
  []

  # Pressure Poisson equation
  [p_diffusion]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
  []

  # Alpha transport: d(alpha)/dt + div(alpha*u) = 0
  # Using rho_mix as coefficient for consistency with the mass flux from RC
  [alpha_1_time]
    type = LinearFVTimeDerivative
    variable = alpha_1
    factor = rho_mix
  []
  [alpha_1_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha_1
    rhie_chow_user_object = 'rc'
    c_alpha = 0.0
    advected_interp_method = ${advected_interp_method}
    limiter_method = ${limiter_method}
    use_nonorthogonal_correction = false
    activate_mules = true
  []
[]

[LinearFVBCs]
  [walls-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom top'
    variable = vel_x
    functor = 0.0
  []
  [walls-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom top'
    variable = vel_y
    functor = 0.0
  []
  [pressure-extrapolation]
    type = LinearFVExtrapolatedPressureBC
    boundary = 'left right bottom top'
    variable = pressure
    use_two_term_expansion = true
  []
  [walls_alpha_1]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha_1
    functor = 0.0
    boundary = 'left right bottom top'
  []
[]

[ICs]
  [alpha_1]
    type = FunctionIC
    variable = 'alpha_1'
    function = alpha_1_init
  []
[]

[Functions]
  [alpha_1_init]
    type = ParsedFunction
    expression = '0.5 * (1.0 - tanh((sqrt((x-${xc})^2 + (y-${yc})^2) - ${R}) / ${eps}))'
  []
[]

[Executioner]
  type = PIMPLEMultiPhase

  number_of_phases = 1

  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  momentum_l_tol = 1e-12
  pressure_l_tol = 1e-12

  rhie_chow_user_objects = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  phase_systems = 'alpha_1_system'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  phase_equation_relaxation = 0.9

  num_iterations = 100

  pressure_absolute_tolerance = 1e-11
  momentum_absolute_tolerance = 1e-11
  phase_absolute_tolerance = 1e-11

  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  phase_petsc_options_iname = '-pc_type -pc_hypre_type'
  phase_petsc_options_value = 'hypre boomeramg'

  print_fields = false
  continue_on_max_its = true
  dt = 1e-3
  num_steps = 2000
  num_piso_iterations = 0

  pin_pressure = true
  pressure_pin_value = 0.0
  pressure_pin_point = '0 0 0'

  enforce_phase_sum = false
  activate_interface_shapening = true
  shapening_type = 'heaviside'
  smoothing_constant = 100.0
[]

[Postprocessors]
  [area_alpha_1]
    type = ElementIntegralVariablePostprocessor
    variable = alpha_1
  []
  [p_inside]
    type = ElementAverageValue
    variable = pressure
    block = 1
  []
  [p_outside]
    type = ElementAverageValue
    variable = pressure
    block = 2
  []
  [delta_p_meas]
    type = ParsedPostprocessor
    expression = 'p_inside - p_outside'
    pp_names = 'p_inside p_outside'
  []
  [R_late]
    type = ParsedPostprocessor
    expression = 'sqrt(area_alpha_1 / 3.141592653589793)'
    pp_names = 'area_alpha_1'
  []
  [delta_p_th]
    type = ParsedPostprocessor
    expression = '${surface_tension} / R_late'
    pp_names = 'R_late'
  []
  [laplace_rel_error]
    type = ParsedPostprocessor
    expression = 'abs(delta_p_meas - delta_p_th) / abs(delta_p_th)'
    pp_names = 'delta_p_meas delta_p_th'
  []
  [max_u]
    type = ElementExtremeValue
    variable = vel_x
    value_type = max
  []
  [min_u]
    type = ElementExtremeValue
    variable = vel_x
    value_type = min
  []
  [max_v]
    type = ElementExtremeValue
    variable = vel_y
    value_type = max
  []
  [min_v]
    type = ElementExtremeValue
    variable = vel_y
    value_type = min
  []
  [max_abs_velocity]
    type = ParsedPostprocessor
    expression = 'sqrt(max(abs(max_u), abs(min_u))^2 + max(abs(max_v), abs(min_v))^2)'
    pp_names = 'max_u min_u max_v min_v'
  []
  [dp_change]
    type = ChangeOverTimePostprocessor
    postprocessor = delta_p_meas
    take_absolute_value = true
    execute_on = 'TIMESTEP_END'
  []
[]

[UserObjects]
  [converged]
    type = Terminator
    expression = 'dp_change < 1e-3 * delta_p_th'
    message = 'Laplace pressure converged to steady state'
    error_level = 'WARNING'
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  [out]
    type = Exodus
  []
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]
