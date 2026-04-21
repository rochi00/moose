# Single-phase VOF dam break: Martin & Moyce (1952)
# One velocity field, mixture properties, MULES advection for alpha.
# This is the interFoam-style formulation:
#   - d(rho_mix)/dt + div(rho_mix*u) = 0
#   - d(rho_mix*u)/dt + div(rho_mix*u*u) = -grad(p) + div(mu_mix*grad(u)) + rho_mix*g
#   - d(alpha)/dt + div(alpha*u) + div(alpha*(1-alpha)*u_c) = 0
#
# Consistent mass-momentum transport enabled via rho_1/rho_2 on the RC object.

rho1 = 998.19
rho2 = 1.185
mu_1 = 1e-3
mu_2 = 1.48e-5
gravity = 9.81
timestep = 0.0001

initial_length = 0.05715
domain_dims_x = ${fparse 5*initial_length}
domain_dims_y = ${fparse 1.25*initial_length}
dam_x = ${initial_length}
dam_y = ${initial_length}

c_alpha = 0.01
advected_interp_method = 'upwind'
limiter_method = 'vanLeer'
MULES_iterations = 3

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${domain_dims_x}'
    dy = '${domain_dims_y}'
    ix = '200'
    iy = '50'
  []
  uniform_refine = 0
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
    # Consistent mass-momentum transport (OpenFOAM interFoam approach)
    rho_1 = ${rho1}
    rho_2 = ${rho2}
    vof_alpha = 'alpha'
    # p_rgh formulation: pressure variable is p_rgh = p - rho*g*h
    gravity_vector = '0 ${fparse -gravity} 0'
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
  # Momentum: d(rho_mix*u)/dt + div(rho_mix*u*u) = -grad(p) + div(mu_mix*grad(u)) + rho_mix*g
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

  # Pressure Poisson
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

  # Alpha transport: d(alpha)/dt + div(phi*alpha) = 0 (purely volumetric, no rho)
  [alpha_time]
    type = LinearFVTimeDerivative
    variable = alpha
  []
  [alpha_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha
    rhie_chow_user_object = 'rc'
    c_alpha = ${c_alpha}
    rho = rho_mix
    advected_interp_method = ${advected_interp_method}
    limiter_method = ${limiter_method}
    use_nonorthogonal_correction = false
    activate_mules = true
    MULES_iterations = ${MULES_iterations}
    # Volumetric alpha equation (uses phi not rho*phi)
    use_volumetric_flux = true
    # Semi-implicit MULES (removes CFL restriction on alpha)
    semi_implicit_mules = true
  []
[]

[LinearFVBCs]
  # Velocity: no-slip on walls, free outflow at top
  [walls_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom'
    variable = vel_x
    functor = 0.0
  []
  [walls_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom'
    variable = vel_y
    functor = 0.0
  []
  [outlet_u]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x
    use_two_term_expansion = false
    boundary = 'top'
  []
  [outlet_v]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y
    use_two_term_expansion = false
    boundary = 'top'
  []

  # Pressure: atmospheric at top, extrapolated on walls
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top'
    variable = pressure
    functor = 0.0
  []
  [pressure_extrapolation]
    type = LinearFVExtrapolatedPressureBC
    boundary = 'left right bottom'
    variable = pressure
    use_two_term_expansion = true
  []

  # Alpha: zero flux on walls, outflow at top
  [walls_alpha]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha
    functor = 0.0
    boundary = 'left right bottom'
  []
  [outlet_alpha]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = alpha
    use_two_term_expansion = false
    boundary = 'top'
  []
[]

[ICs]
  [alpha_ic]
    type = FunctionIC
    variable = 'alpha'
    function = alpha_init
  []
  [pressure_ic]
    type = FunctionIC
    variable = 'pressure'
    function = p_rgh_init
  []
[]

[Functions]
  # Sharp alpha IC
  [alpha_init]
    type = ParsedFunction
    expression = 'if((x < ${dam_x} & y < ${dam_y}), 1.0, 0.0)'
  []
  # p_rgh initial condition: p_rgh = p + rho*g*y (note: g points down, gh = -g*y)
  # p_hydrostatic = rho_mix * g * max(H - y, 0) where H = dam_y within the column
  # p_rgh = p_hydro + rho_mix * g * y = rho_mix * g * H (constant within each phase region)
  # In water column: p_rgh = rho_w*g*dam_y + rho_a*g*(domain_y - dam_y)
  # In air (above or beside): p_rgh = rho_a * g * domain_y (above) or rho_a*g*domain_y (beside)
  # Simplified: p_rgh = alpha * rho_w * g * dam_y + (1-alpha) * rho_a * g * domain_y + blend
  # Exact piecewise for dam geometry:
  [p_rgh_init]
    type = ParsedFunction
    expression = 'a := if((x < ${dam_x} & y < ${dam_y}), 1.0, 0.0);
                  p_rgh := a * (${rho1} * ${gravity} * ${dam_y} + ${rho2} * ${gravity} * (${domain_dims_y} - ${dam_y}))
                         + (1.0 - a) * ${rho2} * ${gravity} * ${domain_dims_y};
                  p_rgh'
  []
[]

[Executioner]
  type = PIMPLEMultiPhase

  number_of_phases = 1

  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  phase_l_abs_tol = 1e-12
  momentum_l_tol = 1e-12
  pressure_l_tol = 1e-12
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
  dt = ${timestep}
  end_time = 0.4
  num_piso_iterations = 0

  enforce_phase_sum = false
  activate_interface_shapening = true
  shapening_type = 'heaviside'
  smoothing_constant = 100.0
  MULES_iterations = ${MULES_iterations}

  # Priority 3: Alpha sub-cycling (allows larger momentum dt)
  n_alpha_subcycles = 2
[]

[AuxVariables]
  [water_heights]
    type = MooseLinearVariableFVReal
  []
  [water_lengths]
    type = MooseLinearVariableFVReal
  []
[]

[AuxKernels]
  [compute_water_heights]
    type = ParsedAux
    variable = 'water_heights'
    coupled_variables = 'alpha'
    expression = 'if(alpha>0.5,y,0)'
    use_xyzt = true
    execute_on = 'TIMESTEP_END'
  []
  [compute_water_lengths]
    type = ParsedAux
    variable = 'water_lengths'
    coupled_variables = 'alpha'
    expression = 'if(alpha>0.5,x,0)'
    use_xyzt = true
    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  [front_height]
    type = ElementExtremeValue
    variable = 'water_heights'
    execute_on = 'TIMESTEP_END'
  []
  [front_length]
    type = ElementExtremeValue
    variable = 'water_lengths'
    execute_on = 'TIMESTEP_END'
  []
  [total_alpha]
    type = ElementIntegralVariablePostprocessor
    variable = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
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
  [max_u]
    type = ElementExtremeValue
    variable = vel_x
    value_type = max
  []
  [max_v]
    type = ElementExtremeValue
    variable = vel_y
    value_type = max
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
