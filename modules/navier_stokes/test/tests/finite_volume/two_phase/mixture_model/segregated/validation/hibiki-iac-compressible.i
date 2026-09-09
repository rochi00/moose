##########################################################
# Vertical upward air-water bubbly flow in the 50.8 mm pipe of Hibiki, Ishii and Xiao (2001),
# for validation of the drift-flux model against the drift velocity measurements compiled in
# Hibiki and Ishii, Int. J. Heat Mass Transfer 45 (2002) 707-721.
#
# Two-dimensional axisymmetric: x is the axial coordinate and the axis of symmetry, y is the
# radius. Flow is upward, so gravity points along -x. The three axial stations of the experiment are
# z/D = 6.00, 30.3 and 53.5. With 120 axial cells over 60 diameters the cell width is exactly D/2,
# so every station at a multiple of half a diameter lands on a face: 6.00 and 53.5 are cut exactly.
# 30.3 is not, and cannot be without a mesh five times finer, since putting all three on faces needs
# a cell width of D/10. It is bracketed instead by sidesets at 29.5, 30.0 and 30.5, and reported by
# linear interpolation between 30.0 and 30.5; the third point bounds the curvature that
# interpolation discards. A sideset at 54.0 is kept as well, so that the bias an earlier version of
# this input carried by reporting the outlet half a diameter downstream of the probe is measured
# rather than estimated.
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
#
# The settings that reproduce the published table live in run_iac_sweep.sh alongside this file, not
# in the defaults here: the comparison is a sweep over three dispersed phase density models and four
# void fractions, so no single set of defaults covers it. The script supplies the distorted particle
# drag law with a swarm exponent of 3/4, a bubble diameter of 2.5 mm, two thousand steps, and the
# relative linear tolerances, and records why each is needed. Run it rather than this input alone.
##########################################################

# Air and water at 0.1 MPa
rho_f = 998.0
mu_f = 1.0e-3
# rho_g is the gas density at the reference pressure p_ref, which is the pressure the outlet
# boundary condition holds. The density at any other point follows from the solved pressure.
rho_g = 1.2
p_ref = 101325
mu_g = 1.8e-5
# Weight of the mixture column assumed by the 'imposed' density model below. Unused by the others.
rho_m_col = 948

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
# Turbulent dissipation rate per unit mass, Hibiki and Ishii (2002) Eq. (20),
#   <eps> = <j>/rho_m * (-dP/dz)_F
# with the frictional gradient from a Blasius factor at the liquid superficial velocity. Prescribed
# rather than solved because this case carries no turbulence model.
eps_diss = 0.0326
# Interfacial area concentration imposed at the inlet. Set to the value measured at z/D = 6.00 so
# that what is tested is the axial development the transport equation predicts between the
# experimental stations, not a guess at the inlet bubble size.
ai_in = 93.2
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
  # Cross-sections on which the area averages are taken. zD6 and zD53p5 are the measured stations;
  # zD29p5, zD30 and zD30p5 bracket the measured station at 30.3; zD54 is the legacy outlet station,
  # retained only to quantify the offset it introduced.
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
  [station_29p5]
    type = ParsedGenerateSideset
    input = station_30
    combinatorial_geometry = 'abs(x - ${fparse 29.5 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD29p5'
  []
  [station_30p5]
    type = ParsedGenerateSideset
    input = station_29p5
    combinatorial_geometry = 'abs(x - ${fparse 30.5 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD30p5'
  []
  [station_53p5]
    type = ParsedGenerateSideset
    input = station_30p5
    combinatorial_geometry = 'abs(x - ${fparse 53.5 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD53p5'
  []
  [station_54]
    type = ParsedGenerateSideset
    input = station_53p5
    combinatorial_geometry = 'abs(x - ${fparse 54 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD54'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system phi_system xi_system'
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

      phase_2_density_name = 'rho_g_var'
      phase_2_density_time_derivative = 'drho_g_dt'
      # d(rho_g)/dp for the isothermal ideal gas below. Supplying it is what lets the pressure
      # driven part of d(rho_m)/dt be assembled on the matrix rather than lagged; without it the
      # outer iteration loses contraction at these void fractions.
      phase_2_density_pressure_derivative = '${fparse rho_g / p_ref}'
      phase_2_viscosity_name = ${mu_g}
      phase_2_specific_heat_name = ${cp_g}
      phase_2_thermal_conductivity_name = ${k_g}

      # The gas density follows the solved pressure, so rho_m varies in time and mixture continuity
      # needs its storage term. Without it the pressure equation imposes div(rho_m u_m) = 0 at every
      # instant, which is a divergence free mass flux rather than a mass balance: the expansion the
      # bubbles undergo as they rise has nowhere to go until the transient has died out. The
      # converged answer is unchanged, since d(rho_m)/dt vanishes there, but only with this on is
      # the march to it a mass balance.
      add_mass_density_transient = true

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
[Variables]
  [interface_area]
    type = MooseLinearVariableFVReal
    solver_sys = xi_system
    initial_condition = ${ai_in}
  []
[]

[LinearFVBCs]
  # The interfacial area equation carries its own boundary conditions; the Physics only generates
  # them for the variables it owns. The inlet value is 6 alpha / d_p for spherical bubbles.
  [xi_inlet]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = interface_area
    boundary = 'left'
    functor = ${ai_in}
  []
  [xi_outlet]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = interface_area
    boundary = 'right'
  []
  [xi_wall]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = interface_area
    boundary = 'top bottom'
  []
[]

[LinearFVKernels]
  # Interfacial area concentration, Fluent 16.4-27 in conservative form
  [xi_time]
    type = LinearFVTimeDerivative
    variable = interface_area
    factor = 'rho_g_var'
    conservative_form = true
  []
  [xi_advection]
    type = LinearFVScalarAdvection
    variable = interface_area
    advected_interp_method_name = upwind
    rhie_chow_user_object = ins_rhie_chow_interpolator
    density = 'rho_g_var'
    u_slip = vel_drift_x
    v_slip = vel_drift_y
    slip_boundaries = 'left right'
  []
  [xi_source_sink]
    type = LinearWCNSFV2PInterfaceAreaSourceSink
    variable = interface_area
    model = hibiki-ishii
    u = vel_x
    v = vel_y
    rho_d = 'rho_g_var'
    rho_f = ${rho_f}
    mu_f = ${mu_f}
    fraction_dispersed = phase_2
    sigma = ${sigma}
    epsilon = ${eps_diss}
    gravity = '-9.81 0 0'
  []

  [radial_rz_viscous]
    type = LinearFVRZViscousSource
    variable = vel_y
    mu = 'mu_mixture'
    momentum_component = 'y'
  []
[]

[FunctorMaterials]
  # Only one of the three density models below is active. They differ solely in where the pressure
  # that sets the gas density comes from, which is the comparison section 8 of the report makes.
  # Select one from the command line, e.g.
  #   FunctorMaterials/active='rho_g_imposed c_d j_axial alpha_j alpha_v_gj'
  active = 'rho_g_solved p_dot drho_g_dt_solved c_d j_axial alpha_j alpha_v_gj'

  # Gas density from the solved pressure, isothermal ideal gas: rho_g = rho_g(p_ref) p_abs / p_ref.
  #
  # 'pressure' is the total pressure the flow Physics solves, not a dynamic pressure: gravity enters
  # the momentum equation as a full rho_m g body force (solve_for_dynamic_pressure defaults to
  # false), so the solved field already carries the hydrostatic head, gauged to zero at the outlet.
  # The absolute pressure is therefore p_ref + pressure.
  #
  # This closes the loop that an imposed rho_g(x) leaves open. The column is weighed by whatever
  # rho_mixture the solution produces rather than by an assumed constant, so the expansion is
  # computed rather than prescribed. The flow solver is still incompressible -- the volume created
  # by the expansion does not accelerate u_m -- but the conservative phase equation carries rho_d,
  # so dispersed phase mass is conserved and alpha grows as rho_g falls.
  #
  # The floor guards the first outer iterations, where the pressure field is still developing from
  # its zero initial condition and can undershoot; it is inactive in the converged solution.
  [rho_g_solved]
    type = ParsedFunctorMaterial
    property_name = 'rho_g_var'
    expression = '${rho_g} * max(${p_ref} + pressure, 0.1 * ${p_ref}) / ${p_ref}'
    functor_names = 'pressure'
  []

  # d(rho_d)/dt for the phase equation. rho_d depends on time only through the pressure, so the
  # chain rule gives d(rho_d)/dt = (d rho_d/dp) dp/dt, and dp/dt comes from the pressure variable's
  # own time derivative, which is built from the time integrator rather than differenced by hand.
  [p_dot]
    type = GenericFunctorTimeDerivativeMaterial
    prop_names = 'p_dot'
    prop_values = 'pressure'
  []
  [drho_g_dt_solved]
    type = ParsedFunctorMaterial
    property_name = 'drho_g_dt'
    expression = '${rho_g} / ${p_ref} * p_dot'
    functor_names = 'p_dot'
  []
  # The imposed and constant density models carry no time dependence at all: 'imposed' is a
  # function of height and 'constant' is a number, so d(rho_d)/dt vanishes for both. The density
  # time derivative has to be selected together with the density itself, exactly as every other
  # term of the phase and interfacial area equations does.
  [drho_g_dt_zero]
    type = ParsedFunctorMaterial
    property_name = 'drho_g_dt'
    expression = '0'
  []

  # Gas density from a prescribed hydrostatic column of fixed weight. This presupposes the void
  # fraction, since rho_m_col is the mixture density the answer is supposed to produce, and it gives
  # every void fraction the same expansion ratio. Retained for comparison only.
  [rho_g_imposed]
    type = ParsedFunctorMaterial
    property_name = 'rho_g_var'
    expression = '${rho_g} * (${p_ref} + ${rho_m_col} * 9.81 * (${L} - x)) / ${p_ref}'
  []

  # Gas density held at its reference value. The bubbles do not expand, which isolates how much of
  # the measured interfacial area development the coalescence and breakage terms can supply alone.
  [rho_g_constant]
    type = ParsedFunctorMaterial
    property_name = 'rho_g_var'
    expression = '${rho_g}'
  []

  # The drift-flux quantities of interest, built up in order. Each is a plain functor of the
  # solved fields and of the slip velocity the closure produces.
  [c_d]
    type = ParsedFunctorMaterial
    property_name = 'c_d'
    expression = 'phase_2 * rho_g_var / rho_mixture'
    functor_names = 'phase_2 rho_g_var rho_mixture'
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
  [alpha_54]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD54'
  []
  [j_54]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD54'
  []
  [alpha_j_54]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD54'
  []
  [alpha_v_gj_54]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD54'
  []
  # C0 = <alpha j> / (<alpha><j>), Vgj = <alpha (1-alpha) u_slip> / <alpha>
  [C0_54]
    type = ParsedPostprocessor
    expression = 'alpha_j_54 / (alpha_54 * j_54)'
    pp_names = 'alpha_j_54 alpha_54 j_54'
  []
  [Vgj_54]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_54 / alpha_54'
    pp_names = 'alpha_v_gj_54 alpha_54'
  []
  [alpha_30p5]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD30p5'
  []
  [j_30p5]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD30p5'
  []
  [alpha_j_30p5]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD30p5'
  []
  [alpha_v_gj_30p5]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD30p5'
  []
  # C0 = <alpha j> / (<alpha><j>), Vgj = <alpha (1-alpha) u_slip> / <alpha>
  [C0_30p5]
    type = ParsedPostprocessor
    expression = 'alpha_j_30p5 / (alpha_30p5 * j_30p5)'
    pp_names = 'alpha_j_30p5 alpha_30p5 j_30p5'
  []
  [Vgj_30p5]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_30p5 / alpha_30p5'
    pp_names = 'alpha_v_gj_30p5 alpha_30p5'
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
  # Interfacial area at the bracketing station used only to bound the curvature that the
  # linear interpolation to z/D = 30.3 discards.
  [ai_29p5]
    type = SideAverageValue
    variable = interface_area
    boundary = 'zD29p5'
  []
  [ai_30p5]
    type = SideAverageValue
    variable = interface_area
    boundary = 'zD30p5'
  []
  [ai_53p5]
    type = SideAverageValue
    variable = interface_area
    boundary = 'zD53p5'
  []
[]

[Executioner]
  type = PIMPLE
  num_steps = 400
  dt = 0.02
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'

  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  active_scalar_systems = 'phi_system xi_system'

  momentum_equation_relaxation = 0.7
  active_scalar_equation_relaxation = '0.7 0.7'
  pressure_variable_relaxation = 0.3

  num_iterations = 30
  pressure_absolute_tolerance = 1e-9
  momentum_absolute_tolerance = 1e-9
  active_scalar_absolute_tolerance = '1e-9 1e-9'
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

[Postprocessors]
  [ai_6]
    type = SideAverageValue
    variable = interface_area
    boundary = 'zD6'
  []
  [ai_30]
    type = SideAverageValue
    variable = interface_area
    boundary = 'zD30'
  []
  [ai_54]
    type = SideAverageValue
    variable = interface_area
    boundary = 'zD54'
  []
  [dsm_6]
    type = ParsedPostprocessor
    expression = '6 * alpha_6 / ai_6'
    pp_names = 'alpha_6 ai_6'
  []
  [dsm_30]
    type = ParsedPostprocessor
    expression = '6 * alpha_30 / ai_30'
    pp_names = 'alpha_30 ai_30'
  []
  [dsm_54]
    type = ParsedPostprocessor
    expression = '6 * alpha_54 / ai_54'
    pp_names = 'alpha_54 ai_54'
  []

  # Evidence that the column is weighed by the solution rather than by an assumed constant. The
  # inlet gauge pressure is the weight of the mixture above it, so it falls as the void fraction
  # rises; rho_g_ratio is the expansion the gas undergoes between the inlet and the outlet, which
  # is the quantity that drives the void and interfacial area development of section 8.
  [p_inlet]
    type = SideAverageValue
    variable = pressure
    boundary = 'left'
  []
  [rho_g_inlet]
    type = SideAverageFunctorPostprocessor
    functor = rho_g_var
    boundary = 'left'
  []
  [rho_g_outlet]
    type = SideAverageFunctorPostprocessor
    functor = rho_g_var
    boundary = 'right'
  []
  [rho_g_ratio]
    type = ParsedPostprocessor
    expression = 'rho_g_inlet / rho_g_outlet'
    pp_names = 'rho_g_inlet rho_g_outlet'
  []
[]
