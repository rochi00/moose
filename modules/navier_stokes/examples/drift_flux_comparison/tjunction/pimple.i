# The linear finite volume drift-flux model run as a transient, see simple.i for the model
!include common.i

[Problem]
  linear_sys_names = 'u_system v_system w_system pressure_system phi_system'
[]

[Physics]
  [NavierStokes]
    [FlowSegregated]
      [flow]
        velocity_variable = 'vel_x vel_y vel_z'
        pressure_variable = pressure

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
    [TwoPhaseMixtureSegregated]
      [mixture]
        system_names = phi_system
        phase_1_fraction_name = phase_1
        phase_2_fraction_name = phase_2

        # Phase transport equation, advected at the drift velocity
        add_phase_transport_equation = true
        phase_advection_interpolation = upwind
        phase_fraction_inlet_type = fixed-value
        phase_fraction_inlet_functors = ${inlet_phase_2}

        # Drift flux: the diffusion stress in the momentum equations, with the slip velocity
        # closed by the drag model solved together with the force balance
        add_drift_flux_momentum_terms = true
        density_interp_method = average
        use_dispersed_phase_drag_model = true
        particle_diameter = ${dp}
        add_advection_slip_term = false

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

[Executioner]
  # Same time march as newton.i, so that the two are compared on equal terms
  type = PIMPLE
  dt = 10
  end_time = 4000
  steady_state_detection = true
  steady_state_tolerance = 1e-8
  momentum_systems = 'u_system v_system w_system'
  pressure_system = pressure_system
  active_scalar_systems = phi_system
  rhie_chow_user_object = ins_rhie_chow_interpolator

  momentum_equation_relaxation = 0.5
  pressure_variable_relaxation = 0.2
  active_scalar_equation_relaxation = 0.7
  num_iterations = 50
  continue_on_max_its = true
  momentum_absolute_tolerance = 1e-8
  pressure_absolute_tolerance = 1e-8
  active_scalar_absolute_tolerance = 1e-8

  momentum_l_tol = 0
  pressure_l_tol = 0
  active_scalar_l_tol = 0
  momentum_l_abs_tol = 1e-10
  pressure_l_abs_tol = 1e-10
  active_scalar_l_abs_tol = 1e-10
  momentum_petsc_options_iname = '-pc_type'
  momentum_petsc_options_value = lu
  pressure_petsc_options_iname = '-pc_type'
  pressure_petsc_options_value = lu
  active_scalar_petsc_options_iname = '-pc_type'
  active_scalar_petsc_options_value = lu
[]
