##########################################################
# Vertical upward air-water bubbly flow in the 50.8 mm pipe of Hibiki, Ishii and Xiao (2001),
# for validation of the drift-flux model against the drift velocity measurements compiled in
# Hibiki and Ishii, Int. J. Heat Mass Transfer 45 (2002) 707-721.
#
# Two-dimensional axisymmetric: x is the axial coordinate and the axis of symmetry, y is the
# radius. Flow is upward, so gravity points along -x. The three axial stations of the experiment,
# z/D = 6.00, 30 and 54, are cut as internal sidesets; with 120 axial cells over 60 diameters they
# land exactly on cell faces.
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

# Geometry
D = 0.0508
R = 0.0254
L = 3.048       # 60 diameters

# Flow condition, the first of the low void fraction cases of the reference report
alpha_in = 0.049
jf = 0.491
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
    nx = 120
    ny = 16
  []
  # Cross-sections at z/D = 6, 30 and 54, on which the area averages are taken
  [station_6]
    type = ParsedGenerateSideset
    input = gen
    combinatorial_geometry = 'abs(x - ${fparse 6 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD6'
  []
  [station_30]
    type = ParsedGenerateSideset
    input = station_6
    combinatorial_geometry = 'abs(x - ${fparse 30 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD30'
  []
  [station_53p5]
    type = ParsedGenerateSideset
    input = station_30
    combinatorial_geometry = 'abs(x - ${fparse 53.5 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD53p5'
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
  [alpha_6]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD6'
  []
  [j_6]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD6'
  []
  [alpha_j_6]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD6'
  []
  [alpha_v_gj_6]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD6'
  []
  # C0 = <alpha j> / (<alpha><j>), Vgj = <alpha (1-alpha) u_slip> / <alpha>
  [C0_6]
    type = ParsedPostprocessor
    expression = 'alpha_j_6 / (alpha_6 * j_6)'
    pp_names = 'alpha_j_6 alpha_6 j_6'
  []
  [Vgj_6]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_6 / alpha_6'
    pp_names = 'alpha_v_gj_6 alpha_6'
  []
  [alpha_30]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD30'
  []
  [j_30]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD30'
  []
  [alpha_j_30]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD30'
  []
  [alpha_v_gj_30]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD30'
  []
  # C0 = <alpha j> / (<alpha><j>), Vgj = <alpha (1-alpha) u_slip> / <alpha>
  [C0_30]
    type = ParsedPostprocessor
    expression = 'alpha_j_30 / (alpha_30 * j_30)'
    pp_names = 'alpha_j_30 alpha_30 j_30'
  []
  [Vgj_30]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_30 / alpha_30'
    pp_names = 'alpha_v_gj_30 alpha_30'
  []
  [alpha_53p5]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD53p5'
  []
  [j_53p5]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD53p5'
  []
  [alpha_j_53p5]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD53p5'
  []
  [alpha_v_gj_53p5]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD53p5'
  []
  # C0 = <alpha j> / (<alpha><j>), Vgj = <alpha (1-alpha) u_slip> / <alpha>
  [C0_53p5]
    type = ParsedPostprocessor
    expression = 'alpha_j_53p5 / (alpha_53p5 * j_53p5)'
    pp_names = 'alpha_j_53p5 alpha_53p5 j_53p5'
  []
  [Vgj_53p5]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_53p5 / alpha_53p5'
    pp_names = 'alpha_v_gj_53p5 alpha_53p5'
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

# Radial profiles at the two experimental stations. This is the quantity Swiderski et al. validate
# against for this same facility, and the one a flat void profile cannot reproduce.
[Positions]
  # The 16 radial cell centroids at each experimental station
  [radial_zD6]
    type = InputPositions
    positions = '0.304800 0.000794 0' '0.304800 0.002381 0' '0.304800 0.003969 0' '0.304800 0.005556 0' '0.304800 0.007144 0' '0.304800 0.008731 0' '0.304800 0.010319 0' '0.304800 0.011906 0' '0.304800 0.013494 0' '0.304800 0.015081 0' '0.304800 0.016669 0' '0.304800 0.018256 0' '0.304800 0.019844 0' '0.304800 0.021431 0' '0.304800 0.023019 0' '0.304800 0.024606 0'
  []
  [radial_zD53p5]
    type = InputPositions
    positions = '2.717800 0.000794 0' '2.717800 0.002381 0' '2.717800 0.003969 0' '2.717800 0.005556 0' '2.717800 0.007144 0' '2.717800 0.008731 0' '2.717800 0.010319 0' '2.717800 0.011906 0' '2.717800 0.013494 0' '2.717800 0.015081 0' '2.717800 0.016669 0' '2.717800 0.018256 0' '2.717800 0.019844 0' '2.717800 0.021431 0' '2.717800 0.023019 0' '2.717800 0.024606 0'
  []
[]

[VectorPostprocessors]
  # Radial void fraction and axial velocity, the quantity Swiderski et al. validate against for
  # this facility. Sampled as functors because the variables live in a linear FV system.
  [prof_zD6]
    type = PositionsFunctorValueSampler
    positions = radial_zD6
    functors = 'phase_2 vel_x'
    discontinuous = true
    sort_by = y
    execute_on = 'FINAL'
  []
  [prof_zD53p5]
    type = PositionsFunctorValueSampler
    positions = radial_zD53p5
    functors = 'phase_2 vel_x'
    discontinuous = true
    sort_by = y
    execute_on = 'FINAL'
  []
[]
