##########################################################
# The 50.8 mm pipe of hibiki-50mm-pipe.i, built instead with the nonlinear (Newton) finite volume
# mixture model. Same mesh, stations, closure (Schiller-Naumann drag solved for the slip speed,
# mixture density in the closure, diffusion velocity in the phase equation) and time march, so that
# the two discretizations can be compared on equal terms, both in the drift-flux parameters they
# produce and in the cost of producing them.
##########################################################

rho_f = 998.0
mu_f = 1.0e-3
rho_g = 1.2
mu_g = 1.8e-5

D = 0.0508
R = 0.0254
L = 3.048

alpha_in = 0.049
jf = 0.491
jg = 0.0275
dp = 0.003

rho_m_in = '${fparse alpha_in * rho_g + (1 - alpha_in) * rho_f}'
u_in = '${fparse (rho_g * jg + rho_f * jf) / rho_m_in}'

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
  # Cross-sections on which the area averages are taken
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
  [station_54]
    type = ParsedGenerateSideset
    input = station_30
    combinatorial_geometry = 'abs(x - ${fparse 54 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD54'
  []
[]

[Physics]
  [NavierStokes]
    [Flow/flow]
      compressibility = 'incompressible'
      density = 'rho_mixture'
      dynamic_viscosity = 'mu_mixture'
      gravity = '-9.81 0 0'

      initial_velocity = '${u_in} 0 0'
      initial_pressure = 0

      inlet_boundaries = 'left'
      momentum_inlet_types = 'fixed-velocity'
      momentum_inlet_functors = '${u_in} 0'

      wall_boundaries = 'top bottom'
      momentum_wall_types = 'noslip symmetry'

      outlet_boundaries = 'right'
      momentum_outlet_types = 'fixed-pressure'
      pressure_functors = '0'

      mass_advection_interpolation = 'upwind'
      momentum_advection_interpolation = 'upwind'
      velocity_interpolation = 'rc'
    []
    [TwoPhaseMixture/mixture]
      phase_1_fraction_name = 'phase_1'
      phase_2_fraction_name = 'phase_2'

      add_phase_transport_equation = true
      phase_advection_interpolation = 'upwind'
      phase_fraction_inlet_type = 'fixed-value'
      phase_fraction_inlet_functors = '${alpha_in}'
      ghost_layers = 5

      add_drift_flux_momentum_terms = true
      density_interp_method = 'average'

      phase_1_density_name = ${rho_f}
      phase_1_viscosity_name = ${mu_f}
      phase_1_specific_heat_name = ${cp_f}
      phase_1_thermal_conductivity_name = ${k_f}

      phase_2_density_name = ${rho_g}
      phase_2_viscosity_name = ${mu_g}
      phase_2_specific_heat_name = ${cp_g}
      phase_2_thermal_conductivity_name = ${k_g}

      use_dispersed_phase_drag_model = true
      particle_diameter = ${dp}
    []
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  num_steps = 400
  dt = 0.02
  nl_rel_tol = 1e-8
  nl_max_its = 30
  line_search = 'basic'
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
  []
[]


[FunctorMaterials]
  # The same drift-flux quantities as hibiki-pipe-base.i, in AD form because the nonlinear slip
  # material declares an AD slip velocity
  [c_d]
    type = ADParsedFunctorMaterial
    property_name = 'c_d'
    expression = 'phase_2 * ${rho_g} / rho_mixture'
    functor_names = 'phase_2 rho_mixture'
  []
  [j_axial]
    type = ADParsedFunctorMaterial
    property_name = 'j_axial'
    expression = 'vel_x + (phase_2 - c_d) * vel_slip_x'
    functor_names = 'vel_x phase_2 c_d vel_slip_x'
  []
  [alpha_j]
    type = ADParsedFunctorMaterial
    property_name = 'alpha_j'
    expression = 'phase_2 * j_axial'
    functor_names = 'phase_2 j_axial'
  []
  [alpha_v_gj]
    type = ADParsedFunctorMaterial
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
  # C0 = <alpha j> / (<alpha><j>), Vgj = <alpha (1-alpha) u_slip> / <alpha>
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
  [avg_alpha]
    type = ElementAverageValue
    variable = phase_2
  []
[]

[Outputs]
  csv = true
[]
