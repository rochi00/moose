##########################################################
# The 50.8 mm pipe of hibiki-50mm-pipe.i, built instead with the nonlinear (Newton) finite volume
# mixture model, to compare the numerical behaviour of the two discretizations on the case that
# defeats the linear FV steady solve.
#
# This is a comparison of NUMERICS, not of physics. The nonlinear implementation carries the
# pre-correction forms throughout: the opposite drift stress sign, the dilute coefficient, the slip
# rather than the diffusion velocity in the phase equation, and the continuous phase density in the
# closure. Its answer is therefore not expected to agree with the linear FV one. What is being
# asked is only whether the phase fraction stays bounded.
##########################################################

rho_f = 998.0
mu_f = 1.0e-3
rho_g = 1.2
mu_g = 1.8e-5

R = 0.0254
L = 3.048

alpha_in = 0.20
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
  type = Steady
  solve_type = 'NEWTON'
  nl_rel_tol = 1e-8
  nl_max_its = 50
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
  [avg_alpha]
    type = ElementAverageValue
    variable = phase_2
  []
[]

[Outputs]
  csv = true
[]
