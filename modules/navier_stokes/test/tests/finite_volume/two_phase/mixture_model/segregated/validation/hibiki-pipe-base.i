##########################################################
# Common part of the vertical upward air-water bubbly pipe flow cases of Hibiki and coworkers,
# for validation of the drift-flux model against the drift velocity measurements compiled in
# Hibiki and Ishii, Int. J. Heat Mass Transfer 45 (2002) 707-721.
#
# Two-dimensional axisymmetric: x is the axial coordinate and the axis of symmetry, y is the
# radius. Flow is upward, so gravity points along -x. The including input sets the geometry (D,
# R, L, nx), the liquid superficial velocity jf, and cuts its own measuring stations as internal
# sidesets with the area averages taken on them.
#
# The quantities of interest are the drift-flux parameters of the one-dimensional model,
#
#   C0   = <alpha j> / (<alpha> <j>),        Vgj = <alpha (1-alpha) u_slip> / <alpha>
#
# formed from the mixture velocity and the slip velocity the closure produces, using
#
#   j    = u_m + (alpha - c_d) u_slip,       c_d = alpha rho_g / rho_m
#
# See the validation section of mms-drift-flux.tex for the derivation and for what these can and
# cannot establish. In short: Vgj is a genuine test of the algebraic slip closure, whereas C0 is
# expected to come out near unity because the model carries no lift force, and is reported as a
# diagnostic of that.
##########################################################

# Air and water at 0.1 MPa
rho_f = 998.0
mu_f = 1.0e-3
rho_g = 1.2
mu_g = 1.8e-5

# Flow condition, the first of the low void fraction cases of the reference report
alpha_in = 0.049
jg = 0.0275
# Sauter mean diameter, within the 1.40 to 3.86 mm range reported at these stations
dp = 0.003
sigma = 0.0728
drag_model = 'schiller-naumann'
swarm_exponent = 0.0

# Mass-averaged mixture velocity at the inlet, (rho_g jg + rho_f jf) / rho_m
rho_m_in = '${fparse alpha_in * rho_g + (1 - alpha_in) * rho_f}'
u_in = '${fparse (rho_g * jg + rho_f * jf) / rho_m_in}'

# Required by the mixture Physics but not used: there is no energy equation here
cp_f = 1
k_f = 1
cp_g = 1
k_g = 1

[Mesh]
  coord_type = RZ
  rz_coord_axis = X
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${L}
    ymin = 0
    ymax = ${R}
    nx = ${nx}
    ny = 16
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system phi_system'
  previous_nl_solution_required = true
[]

[Physics]
  [NavierStokes]
    [FlowSegregated/flow]
      compressibility = 'incompressible'
      density = 'rho_mixture'
      dynamic_viscosity = 'mu_mixture'
      gravity = '-9.81 0 0'

      initial_velocity = '${u_in} 0 0'
      initial_pressure = 0

      inlet_boundaries = 'left'
      momentum_inlet_types = 'fixed-velocity'
      momentum_inlet_functors = '${u_in} 0'

      # 'top' is the pipe wall, 'bottom' is the axis of symmetry
      wall_boundaries = 'top bottom'
      momentum_wall_types = 'noslip symmetry'

      outlet_boundaries = 'right'
      momentum_outlet_types = 'fixed-pressure'
      pressure_functors = '0'

      momentum_advection_interpolation = 'upwind'
      orthogonality_correction = false
      pressure_two_term_bc_expansion = true
      momentum_two_term_bc_expansion = true
    []
    [TwoPhaseMixtureSegregated/mixture]
      system_names = 'phi_system'
      phase_1_fraction_name = 'phase_1'
      phase_2_fraction_name = 'phase_2'

      add_phase_transport_equation = true
      phase_advection_interpolation = 'upwind'
      phase_fraction_inlet_type = 'fixed-value'
      phase_fraction_inlet_functors = '${alpha_in}'

      # Drift flux
      add_drift_flux_momentum_terms = true
      density_interp_method = 'average'

      # Continuous phase is water, dispersed phase is air
      phase_1_density_name = ${rho_f}
      phase_1_viscosity_name = ${mu_f}
      phase_1_specific_heat_name = ${cp_f}
      phase_1_thermal_conductivity_name = ${k_f}

      phase_2_density_name = ${rho_g}
      phase_2_viscosity_name = ${mu_g}
      phase_2_specific_heat_name = ${cp_g}
      phase_2_thermal_conductivity_name = ${k_g}

      use_dispersed_phase_drag_model = true
      slip_drag_model = ${drag_model}
      surface_tension = ${sigma}
      slip_swarm_exponent = ${swarm_exponent}
      particle_diameter = ${dp}
      add_advection_slip_term = false
    []
  []
[]

# The radial momentum equation needs the axisymmetric hoop term, which the Physics does not add
[LinearFVKernels]
  [radial_rz_viscous]
    type = LinearFVRZViscousSource
    variable = vel_y
    mu = 'mu_mixture'
    momentum_component = 'y'
  []
[]

[FunctorMaterials]
  # The drift-flux quantities of interest, built up in order. Each is a plain functor of the
  # solved fields and of the slip velocity the closure produces.
  [c_d]
    type = ParsedFunctorMaterial
    property_name = 'c_d'
    expression = 'phase_2 * ${rho_g} / rho_mixture'
    functor_names = 'phase_2 rho_mixture'
  []
  [j_axial]
    type = ParsedFunctorMaterial
    property_name = 'j_axial'
    expression = 'vel_x + (phase_2 - c_d) * vel_slip_x'
    functor_names = 'vel_x phase_2 c_d vel_slip_x'
  []
  [alpha_j]
    type = ParsedFunctorMaterial
    property_name = 'alpha_j'
    expression = 'phase_2 * j_axial'
    functor_names = 'phase_2 j_axial'
  []
  [alpha_v_gj]
    type = ParsedFunctorMaterial
    property_name = 'alpha_v_gj'
    expression = 'phase_2 * (1 - phase_2) * vel_slip_x'
    functor_names = 'phase_2 vel_slip_x'
  []
[]

[Postprocessors]
  [min_alpha]
    type = ElementExtremeValue
    variable = phase_2
    value_type = min
  []
  [max_alpha]
    type = ElementExtremeValue
    variable = phase_2
  []
[]

[Executioner]
  type = PIMPLE
  num_steps = 400
  dt = 0.02
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'

  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  active_scalar_systems = 'phi_system'

  momentum_equation_relaxation = 0.7
  active_scalar_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3

  num_iterations = 30
  pressure_absolute_tolerance = 1e-9
  momentum_absolute_tolerance = 1e-9
  active_scalar_absolute_tolerance = 1e-9
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  active_scalar_petsc_options_iname = '-pc_type -pc_hypre_type'
  active_scalar_petsc_options_value = 'hypre boomeramg'
  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  active_scalar_l_abs_tol = 1e-12
  continue_on_max_its = true
[]

[Outputs]
  csv = true
[]
