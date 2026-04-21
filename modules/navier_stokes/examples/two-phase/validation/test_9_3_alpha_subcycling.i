# Validation Test 9.3: Alpha Sub-Cycling Equivalence
# Pure advection of a liquid column in a uniform velocity field.
# Based on plic_advection_2d.i but using algebraic MULES advection.
#
# Run this file twice:
#   A) As-is:           dt=0.01,  n_alpha_subcycles=1  (baseline)
#   B) With overrides:  dt=0.04,  n_alpha_subcycles=4  (subcycled)
#
# Command for case B:
#   moose_app -i test_9_3_alpha_subcycling.i Executioner/dt=0.04 Executioner/n_alpha_subcycles=4 Executioner/num_steps=5 Outputs/csv/file_base=test_9_3_subcycled
#
# Success criteria:
#   - total_alpha (mass) identical in both cases at same physical times
#   - alpha_min >= 0, alpha_max <= 1 (boundedness)
#   - Interface shape nearly identical

rho1 = 1.0
rho2 = 1.0
mu_1 = 1e-3
mu_2 = 1e-3

advected_interp_method = 'upwind'
limiter_method = 'vanLeer'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 0.5
    nx = 40
    ny = 20
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
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_condition = 1.0
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

  [alpha_time]
    type = LinearFVTimeDerivative
    variable = alpha
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
  [walls-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
    variable = vel_x
    functor = 1.0
  []
  [walls-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
    variable = vel_y
    functor = 0.0
  []
  [inlet-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = vel_x
    functor = 1.0
  []
  [inlet-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = vel_y
    functor = 0.0
  []
  [outlet-u]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x
    use_two_term_expansion = false
    boundary = 'right'
  []
  [outlet-v]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y
    use_two_term_expansion = false
    boundary = 'right'
  []
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'right'
    variable = pressure
    functor = 0.0
  []
  [pressure-extrapolation]
    type = LinearFVExtrapolatedPressureBC
    boundary = 'left top bottom'
    variable = pressure
    use_two_term_expansion = true
  []
  [walls_alpha]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha
    functor = 0.0
    boundary = 'top bottom'
  []
  [inlet_alpha]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = alpha
    functor = 0.0
    boundary = 'left'
  []
  [outlet_alpha]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = alpha
    use_two_term_expansion = false
    boundary = 'right'
  []
[]

[ICs]
  [alpha_ic]
    type = FunctionIC
    variable = alpha
    function = alpha_init
  []
[]

[Functions]
  [alpha_init]
    type = ParsedFunction
    expression = 'if(x < 0.3, 1.0, 0.0)'
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
  phase_systems = 'alpha_system'

  momentum_equation_relaxation = 1.0
  pressure_variable_relaxation = 1.0
  phase_equation_relaxation = 1.0

  num_iterations = 1

  pressure_absolute_tolerance = 1e-11
  momentum_absolute_tolerance = 1e-11
  phase_absolute_tolerance = 1e-11

  momentum_petsc_options_iname = '-pc_type'
  momentum_petsc_options_value = 'jacobi'
  pressure_petsc_options_iname = '-pc_type'
  pressure_petsc_options_value = 'jacobi'
  phase_petsc_options_iname = '-pc_type'
  phase_petsc_options_value = 'jacobi'

  print_fields = false
  continue_on_max_its = true
  dt = 0.01
  num_steps = 20
  num_piso_iterations = 1

  pin_pressure = false
  enforce_phase_sum = false
  activate_interface_shapening = false

  # Sub-cycling parameter (override to 4 for case B)
  n_alpha_subcycles = 1
[]

[Postprocessors]
  [total_alpha]
    type = ElementIntegralVariablePostprocessor
    variable = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [alpha_max]
    type = ElementExtremeValue
    variable = alpha
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [alpha_min]
    type = ElementExtremeValue
    variable = alpha
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]
