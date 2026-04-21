# Validation Test 9.2: Hydrostatic Column — Consistent Mass-Momentum Transport
# Single-phase VOF formulation with high density ratio (1000:1).
# Heavy fluid on bottom, light fluid on top, gravity pointing down.
#
# Expected at equilibrium:
#   - velocity = 0 everywhere (to machine precision with consistent transport)
#   - hydrostatic pressure: p_bottom ~ rho1*g*h1 + rho2*g*h2
#   - alpha unchanged (sharp, bounded)
#
# Without consistent rhoPhi: spurious velocities ~ O(rho1/rho2)
# With consistent rhoPhi: spurious velocities ~ O(eps_machine)
#
# The key new parameter is rho_1/rho_2 on the RhieChow object,
# which enables consistent mass flux computation.

rho1 = 1000.0
rho2 = 1.0
mu_1 = 1e-3
mu_2 = 1.8e-5
gravity = 9.81

domain_x = 0.1
domain_y = 0.2
interface_y = 0.1
nxy_x = 20
nxy_y = 40

eps = ${fparse 2.0 * domain_y / nxy_y}

advected_interp_method = 'upwind'
limiter_method = 'vanLeer'

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${domain_x}'
    dy = '${domain_y}'
    ix = '${nxy_x}'
    iy = '${nxy_y}'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system alpha_system'
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
    # Enable consistent mass-momentum transport
    rho_1 = ${rho1}
    rho_2 = ${rho2}
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
  [alpha]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_system
  []
[]

[FunctorMaterials]
  [mixture]
    type = WCNSLinearFVMixtureFunctorMaterial
    phase_1_names = '${rho1} ${mu_1}'
    phase_2_names = '${rho2} ${mu_2}'
    prop_names = 'rho_mix mu_mix'
    phase_1_fraction = alpha
    limit_phase_fraction = true
  []
[]

[LinearFVKernels]
  # Momentum
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
  [v_gravity]
    type = LinearFVSource
    variable = vel_y
    source_density = 'rho_mix'
    scaling_factor = ${fparse -gravity}
  []

  # Pressure
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

  # Alpha transport
  [alpha_time]
    type = LinearFVTimeDerivative
    variable = alpha
    factor = rho_mix
  []
  [alpha_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha
    rhie_chow_user_object = 'rc'
    c_alpha = 0.0
    advected_interp_method = ${advected_interp_method}
    limiter_method = ${limiter_method}
    use_nonorthogonal_correction = false
    activate_mules = true
  []
[]

[LinearFVBCs]
  [walls_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom top'
    variable = vel_x
    functor = 0.0
  []
  [walls_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom top'
    variable = vel_y
    functor = 0.0
  []
  [top_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top'
    variable = pressure
    functor = 0.0
  []
  [walls_p]
    type = LinearFVExtrapolatedPressureBC
    boundary = 'left right bottom'
    variable = pressure
    use_two_term_expansion = true
  []
  [walls_alpha]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha
    functor = 0.0
    boundary = 'left right bottom top'
  []
[]

[ICs]
  [alpha_ic]
    type = FunctionIC
    variable = 'alpha'
    function = alpha_init
  []
[]

[Functions]
  [alpha_init]
    type = ParsedFunction
    expression = '0.5 * (1.0 - tanh((y - ${interface_y}) / ${eps}))'
  []
[]

[Executioner]
  type = PIMPLEMultiPhase

  number_of_phases = 1

  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  momentum_l_tol = 1e-12
  pressure_l_tol = 1e-12
  phase_l_abs_tol = 1e-12
  phase_l_tol = 1e-12

  rhie_chow_user_objects = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  phase_systems = 'alpha_system'

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
  dt = 1e-4
  num_steps = 200

  pin_pressure = false
  enforce_phase_sum = false
  activate_interface_shapening = false
[]

[Postprocessors]
  # Mass conservation
  [total_alpha]
    type = ElementIntegralVariablePostprocessor
    variable = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Boundedness
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
  # Spurious velocities (should be ~0)
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
  # Pressure at bottom (expected: rho1*g*0.1 + rho2*g*0.1 ~ 982 Pa)
  [p_bottom]
    type = SideAverageValue
    variable = pressure
    boundary = 'bottom'
    execute_on = 'TIMESTEP_END'
  []
  [p_bottom_expected]
    type = ParsedPostprocessor
    expression = '${rho1} * ${gravity} * ${interface_y} + ${rho2} * ${gravity} * (${domain_y} - ${interface_y})'
    pp_names = ''
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]
