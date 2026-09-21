!include common.i

[Physics]
  [NavierStokes]
    [Flow]
      [flow]
        compressibility = incompressible
        density = rho_mixture
        dynamic_viscosity = mu_mixture

        gravity = ${gravity}
        solve_for_dynamic_pressure = true
        reference_pressure_point = '1 0 0'
        initial_pressure = 0
        initial_velocity = '1e-6 0 0'

        inlet_boundaries = inlet
        momentum_inlet_types = fixed-velocity
        momentum_inlet_functors = '${inlet_velocity} 0 0'

        outlet_boundaries = 'outlet_minus_y outlet_plus_y'
        momentum_outlet_types = 'fixed-pressure-zero-gradient fixed-pressure-zero-gradient'
        pressure_functors = '0 0'

        wall_boundaries = walls
        momentum_wall_types = noslip

        momentum_advection_interpolation = upwind
      []
    []
    [TwoPhaseMixture]
      [mixture]
        phase_1_fraction_name = phase_1
        phase_2_fraction_name = phase_2

        # Phase transport equation, advected at the slip velocity
        add_phase_transport_equation = true
        phase_advection_interpolation = upwind
        phase_fraction_inlet_type = fixed-value
        phase_fraction_inlet_functors = ${inlet_phase_2}
        ghost_layers = 5

        # Drift flux: the diffusion stress in the momentum equations, with the slip velocity
        # closed by the drag model
        add_drift_flux_momentum_terms = true
        density_interp_method = average
        use_dispersed_phase_drag_model = true
        particle_diameter = ${dp}

        phase_1_density_name = ${rho}
        phase_1_viscosity_name = ${mu}
        phase_1_specific_heat_name = ${cp}
        phase_1_thermal_conductivity_name = ${k}

        phase_2_density_name = ${rho_d}
        phase_2_viscosity_name = ${mu_d}
        phase_2_specific_heat_name = ${cp_d}
        phase_2_thermal_conductivity_name = ${k_d}
        output_all_properties = true
      []
    []
  []
[]

[Preconditioning]
  # The same full single-matrix preconditioner NEWTON builds by itself, spelled out so that the
  # factorization can be handed to MUMPS from the command line, see newton-pipe.i
  [SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu NONZERO'
  []
[]

[Executioner]
  # The phase transport equation has no diagonal of its own in a steady solve, and Newton does
  # not converge on it; marched in time it does, and the transient is run to its steady state.
  # One flow-through of the main channel is 200 s.
  type = Transient
  dt = 10
  end_time = 4000
  steady_state_detection = true
  steady_state_tolerance = 1e-8

  solve_type = NEWTON
  line_search = bt
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 50
[]
