##########################################################
# Vertical upward air-water bubbly flow in the annulus of Hibiki, Situ, Mi and Ishii,
# Int. J. Heat Mass Transfer 46(8) 1479-1496 (2003).
#
# Why this case. The distribution parameter C0 = <alpha j>/(<alpha><j>) departs from unity only
# when void and velocity are radially correlated, and that correlation is produced by lift and wall
# lubrication, neither of which exists in this model. It therefore returns C0 -> 1 by construction.
# This paper reports the phase distribution pattern for every condition it ran, and identifies a
# 'bubble-mixing' region at low liquid flux in which
#
#   "the void fraction profiles were almost uniform along the channel radius"   (for <jf> = 0.272
#   m/s and void fraction lower than 0.10)
#
# A uniform void profile gives C0 = 1 exactly, not by correlation. So at these two conditions the
# model's structural weakness is absent from the measurement as well, and the comparison isolates
# the drift velocity.
#
# NOTE the mechanisms differ. The experiment is flat because of strong bubble-induced turbulent
# mixing; this model is flat because it has no lateral force at all. Agreement on the area averaged
# void fraction is therefore not evidence that the phase distribution is captured, and it should
# not be expected to survive at intermediate liquid flux where the wall peak appears.
#
# Geometry. Inner rod 19.1 mm diameter, outer tube 38.1 mm inner diameter, hydraulic equivalent
# diameter D_H = 19.1 mm. Heater rod overall length 2670 mm. Gas measurements at z/D_H = 40.3 and
# 61.7 averaged to z/D_H = 51.0; the gamma densitometer for area averaged void sits at z/D_H = 51.1.
#
# Conditions, Table 1. <jg,N> is reduced to normal conditions (atmospheric, 20 C), so the local
# superficial gas velocity at the station is lower by the compression of the column above it. The
# gas velocities were chosen to give <alpha> = 0.050, 0.10, 0.15, 0.20, 0.25.
#
#   <jf> = 0.272 m/s :  <jg,N> = 0.0313  0.0506  0.0690  0.0888  0.105
#
# The <alpha> values quoted alongside those in the source (0.050, 0.10, 0.15, 0.20, 0.25) are
# DESIGN TARGETS, not measurements: the gas velocities were "roughly determined so as to provide
# the same area-averaged void fractions". Do not compare against them. Hibiki's own correlations
# miss them by up to a third at the higher gas flows, and the measured values differ.
#
# The measurements are in Figs. 15 and 17, digitised here. Fig. 17 reports only the two lowest gas
# flows at this liquid flux, so only those two conditions can be compared:
#
#   <jg,N>    alpha_in   u_in     measured <alpha>   measured Vgj   (Fig.17, open circles)
#   0.0313    0.0494     0.2862   0.0521             0.2018
#   0.0506    0.0797     0.2956   0.0894             0.1867
#
# Fig. 15 gives the measured distribution parameter for this facility, eight points spanning
# <DSm>/D_H = 0.091 to 0.164:  C0 = 0.969 to 1.074, mean 1.025. That is what makes the case usable:
# this model returns C0 -> 1 by construction, so a comparison is only meaningful where the
# measured C0 is also near unity. Note that Eq.(16) of Hibiki and Ishii reads 1.03 to 1.16 over the
# same range, so the correlation is not a substitute for the measurement here; the source itself
# reports a 10.2 per cent deviation for it in an annulus.
#
# Run the second condition with
#   alpha_in=0.0797 u_in=0.2956
##########################################################

# Air and water at 20 C, 0.1 MPa
rho_f = 998.0
mu_f = 1.0e-3
rho_g = 1.2
mu_g = 1.8e-5
sigma = 0.0728
p_ref = 101325

# Annulus
R_rod = 0.00955
R_tube = 0.01905
D_H = 0.0191
L = 2.670
z_station = '${fparse 51.0 * D_H}'

# Condition. Defaults are <jg,N> = 0.0313, the <alpha> = 0.050 point of Table 1. The inlet values
# below are the mixture velocity and void fraction that deliver <jf> = 0.272 m/s and that <jg,N> at
# the inlet pressure; the computed 'jgN_station' postprocessor is what must be checked against
# Table 1, not these.
alpha_in = 0.0494
u_in = 0.2862

# The drag law that reproduces Ishii's bubbly drift velocity correlation,
# Vgj = sqrt(2) (g sigma drho / rho_f^2)^(1/4) (1 - alpha)^1.75, which is Eq.(13) of the source.
# One power of (1 - alpha) is the definition Vgj = (1 - alpha) u_slip; the remaining 3/4 is the
# swarm hindrance.
drag_model = 'distorted-particle'
# The swarm exponent is 0.25, not the 0.75 of Ishii's correlation itself. The closure forms its
# relaxation time from the continuous phase viscosity, so the factor of (1 - alpha) that the
# mixture density puts into the buoyancy factor survives into the Stokes-limit slip. The distorted
# particle balance takes its square root, so (1 - alpha)^0.5 is already present before any
# hindrance is applied. The residual that recovers Ishii's (1 - alpha)^0.75 relative velocity is
# therefore 0.25. Verified against the correlation to better than 0.01% for alpha in [0.02, 0.30].
swarm_exponent = 0.25
dp = 0.003

# Required by the mixture Physics but unused: there is no energy equation here
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
    ymin = ${R_rod}
    ymax = ${R_tube}
    nx = 148
    ny = 16
  []
  [station]
    type = ParsedGenerateSideset
    input = gen
    combinatorial_geometry = 'abs(x - ${z_station}) < 1e-3'
    normal = '1 0 0'
    new_sideset_name = 'station'
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

      # Both radial boundaries are walls: 'bottom' is the inner rod, 'top' the outer tube.
      # There is no symmetry axis in an annulus.
      wall_boundaries = 'top bottom'
      momentum_wall_types = 'noslip noslip'

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

      add_drift_flux_momentum_terms = true
      density_interp_method = 'average'

      phase_1_density_name = ${rho_f}
      phase_1_viscosity_name = ${mu_f}
      phase_1_specific_heat_name = ${cp_f}
      phase_1_thermal_conductivity_name = ${k_f}

      phase_2_density_name = 'rho_g_var'
      phase_2_density_time_derivative = 'drho_g_dt'
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

[LinearFVKernels]
  [radial_rz_viscous]
    type = LinearFVRZViscousSource
    variable = vel_y
    mu = 'mu_mixture'
    momentum_component = 'y'
  []
[]

[FunctorMaterials]
  # Gas density from the solved pressure. 'pressure' is the total pressure: gravity enters as a
  # full rho_m g body force, so the solved field carries the hydrostatic head, gauged to zero at
  # the outlet. The column is weighed by whatever mixture density the solution produces, which is
  # what makes the expansion of <jg,N> to the local <jg> a computed quantity rather than an input.
  [rho_g_solved]
    type = ParsedFunctorMaterial
    property_name = 'rho_g_var'
    expression = '${rho_g} * max(${p_ref} + pressure, 0.1 * ${p_ref}) / ${p_ref}'
    functor_names = 'pressure'
  []
  [p_dot]
    type = GenericFunctorTimeDerivativeMaterial
    prop_names = 'p_dot'
    prop_values = 'pressure'
  []
  [drho_g_dt]
    type = ParsedFunctorMaterial
    property_name = 'drho_g_dt'
    expression = '${rho_g} / ${p_ref} * p_dot'
    functor_names = 'p_dot'
  []

  # Drift-flux quantities of interest, built in order.
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
  [j_gas]
    type = ParsedFunctorMaterial
    property_name = 'j_gas'
    expression = 'phase_2 * (vel_x + (1 - c_d) * vel_slip_x)'
    functor_names = 'phase_2 vel_x c_d vel_slip_x'
  []
  [j_gas_normal]
    type = ParsedFunctorMaterial
    property_name = 'j_gas_normal'
    expression = 'j_gas * rho_g_var / ${rho_g}'
    functor_names = 'j_gas rho_g_var'
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

[Executioner]
  type = PIMPLE
  num_steps = 600
  dt = 0.05
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'

  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  active_scalar_systems = 'phi_system'

  momentum_equation_relaxation = 0.7
  active_scalar_equation_relaxation = '0.7'
  pressure_variable_relaxation = 0.3

  num_iterations = 30
  pressure_absolute_tolerance = 1e-9
  momentum_absolute_tolerance = 1e-9
  active_scalar_absolute_tolerance = '1e-9'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  active_scalar_petsc_options_iname = '-pc_type -pc_hypre_type'
  active_scalar_petsc_options_value = 'hypre boomeramg'
  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  active_scalar_l_abs_tol = 1e-12
  momentum_l_tol = 1e-10
  pressure_l_tol = 1e-10
  active_scalar_l_tol = 1e-10
  continue_on_max_its = true
[]

[Postprocessors]
  # The comparison quantity: area averaged void fraction at z/D_H = 51.0
  [alpha_station]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'station'
  []
  [alpha_in_pp]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'left'
  []
  # jg reduced to normal conditions, which is what Table 1 of the source reports
  [jgN_station]
    type = SideAverageFunctorPostprocessor
    functor = j_gas_normal
    boundary = 'station'
  []
  [jg_station]
    type = SideAverageFunctorPostprocessor
    functor = j_gas
    boundary = 'station'
  []
  [j_station]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'station'
  []
  [alpha_j_station]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'station'
  []
  [alpha_v_gj_station]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'station'
  []
  # C0 = <alpha j>/(<alpha><j>), Vgj = <alpha (1-alpha) u_slip>/<alpha>
  [C0_station]
    type = ParsedPostprocessor
    expression = 'alpha_j_station / (alpha_station * j_station)'
    pp_names = 'alpha_j_station alpha_station j_station'
  []
  [Vgj_station]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_station / alpha_station'
    pp_names = 'alpha_v_gj_station alpha_station'
  []
  [min_alpha]
    type = ElementExtremeValue
    variable = phase_2
    value_type = min
  []
  [max_alpha]
    type = ElementExtremeValue
    variable = phase_2
  []
  [p_station]
    type = SideAverageValue
    variable = pressure
    boundary = 'station'
  []
[]

[Outputs]
  csv = true
[]
