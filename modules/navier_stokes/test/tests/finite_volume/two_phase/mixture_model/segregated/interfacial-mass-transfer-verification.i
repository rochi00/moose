# Verification of the interfacial mass transfer closed on the transported interfacial area.
#
# A closed, uniform, quiescent box of superheated mixture. The mixture temperature starts above
# saturation, so the interfacial mass transfer generates the dispersed phase, absorbing its latent
# heat from the energy equation and driving the temperature back towards saturation. Every field
# stays spatially uniform, so the case reduces to a system of ordinary differential equations with
# an exact invariant to check against.
#
# With a constant dispersed phase density the material derivative of that density vanishes, and
# with both bubble interaction coefficients set to zero the coalescence and breakage sources
# vanish, so the only source of interfacial area is the phase change term. What remains is
#
#   rho_d d(alpha)/dt = Gamma ,        rho_d d(chi)/dt = (2/3) Gamma chi / alpha ,
#
# and dividing one by the other eliminates Gamma entirely:
#
#   d(chi)/d(alpha) = (2/3) chi / alpha   ==>   chi / alpha^(2/3) = constant .
#
# That invariant is the statement that the transfer grows the particles already present at fixed
# number density, which is the assumption the two thirds exponent encodes. It holds whatever the
# transfer rate, the temperature or the properties are, so the postprocessor below is compared
# against its initial value rather than against a computed answer. Getting the exponent, the sign
# of the coupling or the consistency between the two equations wrong all break it.
#
# The transfer rate itself is not prescribed here. It is computed by the Physics from the solved
# interfacial area, and the same computed rate is what the area equation below is given, which is
# the property this case exists to exercise.

# Continuous phase, liquid
rho = 1000.0
mu = 1e-3
cp = 4000.0
k = 0.6

# Dispersed phase. Its density and specific heat are deliberately those of the continuous phase,
# which is what makes this case well posed and what makes it measure one thing.
#
# Equal densities. A closed rigid box cannot accommodate a transfer between phases of different
# density: the mass it holds is fixed, so converting one phase into the other at fixed phase
# densities would demand a change of volume there is nowhere to take. With the two densities equal
# the mixture density is constant, mass and volume are both conserved, and the dilatation the
# relative motion would produce vanishes identically because c_d = alpha.
#
# Equal specific heats. The mixture specific heat is then constant in time, so the energy equation
# needs no time derivative of rho_m cp_m, and the enthalpy the relative motion carries vanishes
# because it is proportional to the difference of the two. What is left driving the energy equation
# is the latent heat of the transfer alone.
#
# The density is constant in time either way, which is what makes the invariant exact rather than
# approximate.
rho_d = 1000.0
mu_d = 1e-5
cp_d = 4000.0
k_d = 0.03

# Saturation state of the transition. The latent heat is far smaller than any real one, chosen so
# that the phase fraction changes by a factor of five over a handful of steps: the invariant is
# checked against a change large enough to distinguish the two thirds exponent from a neighbouring
# one, which a physical latent heat would leave at the level of the solver tolerance.
T_sat = 372.0
h_fg = 1e4

# Initial state. The superheat drives the transfer.
T_initial = 375.0
alpha_initial = 0.05
area_initial = 300.0

t_end = 2.0

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 0.1
    ymin = 0
    ymax = 0.1
    nx = 4
    ny = 4
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system phase_system energy_system xi_system'
[]

[Physics]
  [NavierStokes]
    [FlowSegregated]
      [flow]
        compressibility = 'weakly-compressible'
        density = 'rho_mixture'
        dynamic_viscosity = 'mu_mixture'

        initial_velocity = '1e-12 1e-12 0'
        initial_pressure = 0

        wall_boundaries = 'top left right bottom'
        momentum_wall_types = 'noslip noslip noslip noslip'
        momentum_wall_functors = '0 0; 0 0; 0 0; 0 0'

        orthogonality_correction = false
        pressure_two_term_bc_expansion = true
        momentum_advection_interpolation = 'upwind'
      []
    []
    [FluidHeatTransferSegregated]
      [energy]
        system_names = 'energy_system'
        thermal_conductivity = 'k_mixture'
        specific_heat = 'cp_mixture'
        initial_temperature = ${T_initial}

        energy_wall_boundaries = 'top left right bottom'
        energy_wall_types = 'heatflux heatflux heatflux heatflux'
        energy_wall_functors = '0 0 0 0'
      []
    []
    [TwoPhaseMixtureSegregated]
      [mixture]
        system_names = 'phase_system'
        fluid_heat_transfer_physics = 'energy'

        phase_1_fraction_name = 'phase_1'
        phase_2_fraction_name = 'phase_2'
        initial_phase_fraction = ${alpha_initial}

        phase_1_density_name = ${rho}
        phase_1_viscosity_name = ${mu}
        phase_1_specific_heat_name = ${cp}
        phase_1_thermal_conductivity_name = ${k}

        phase_2_density_name = ${rho_d}
        phase_2_viscosity_name = ${mu_d}
        phase_2_specific_heat_name = ${cp_d}
        phase_2_thermal_conductivity_name = ${k_d}

        # Close the transfer on the solved area instead of prescribing it. The rate this creates is
        # 'mixture_interfacial_mass_transfer_rate', which the area equation below is given.
        interfacial_area = 'interface_area'
        T_saturation = ${T_sat}
        interfacial_latent_heat = ${h_fg}

        use_dispersed_phase_drag_model = true
        particle_diameter = 1e-3
      []
    []
  []
[]

[FunctorMaterials]
  # The factor the energy time derivative is built on. The mixture specific heat the Physics
  # creates is the mass-weighted average, so this product is the mixture enthalpy density per unit
  # temperature.
  [rho_cp]
    type = ParsedFunctorMaterial
    property_name = 'rho_cp'
    expression = 'rho_mixture * cp_mixture'
    functor_names = 'rho_mixture cp_mixture'
  []
[]

[Variables]
  [interface_area]
    type = MooseLinearVariableFVReal
    solver_sys = xi_system
    initial_condition = ${area_initial}
  []
[]

# The interfacial area transport equation, written out here because no Physics assembles it. With
# no flow the advection contributes nothing, so the transient balances the sources.
[LinearFVKernels]
  [area_time]
    type = LinearFVTimeDerivative
    variable = interface_area
    factor = ${rho_d}
    conservative_form = true
  []
  [area_sources]
    type = LinearWCNSFV2PInterfaceAreaSourceSink
    variable = interface_area
    model = 'hibiki-ishii'

    u = vel_x
    v = vel_y

    rho_d = ${rho_d}
    rho_f = ${rho}
    fraction_dispersed = 'phase_2'
    sigma = 0.059
    epsilon = 1e-6

    # The one rate, computed by the Physics from this same variable
    mass_transfer_rate = 'mixture_interfacial_mass_transfer_rate'

    # Zeroing both interaction coefficients leaves the phase change term as the only source, which
    # is what the invariant in the header was derived for. Hibiki and Ishii has no wake
    # entrainment, so this switches off every bubble interaction.
    gamma_c = 0
    gamma_b = 0
  []
[]

[LinearFVBCs]
[]

[Executioner]
  type = PIMPLE
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'

  dt = 0.2
  end_time = ${t_end}

  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  energy_system = 'energy_system'
  active_scalar_systems = 'phase_system xi_system'

  momentum_equation_relaxation = 0.8
  pressure_variable_relaxation = 0.3
  energy_equation_relaxation = 0.9
  active_scalar_equation_relaxation = '0.9 0.9'

  num_iterations = 200
  pressure_absolute_tolerance = 1e-11
  momentum_absolute_tolerance = 1e-11
  energy_absolute_tolerance = 1e-11
  active_scalar_absolute_tolerance = 1e-11

  momentum_l_abs_tol = 1e-13
  pressure_l_abs_tol = 1e-13
  energy_l_abs_tol = 1e-13
  active_scalar_l_abs_tol = 1e-13
  momentum_l_tol = 1e-12
  pressure_l_tol = 1e-12
  energy_l_tol = 1e-12
  active_scalar_l_tol = 1e-12
  continue_on_max_its = true

  pin_pressure = true
  pressure_pin_value = 0.0
  pressure_pin_point = '0.0 0.0 0.0'

  print_fields = false
[]

[Postprocessors]
  [alpha]
    type = ElementAverageValue
    variable = phase_2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [area]
    type = ElementAverageValue
    variable = interface_area
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [temperature]
    type = ElementAverageValue
    variable = T_fluid
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # chi / alpha^(2/3), constant under the phase change term alone. See the header.
  [number_invariant]
    type = ParsedPostprocessor
    expression = 'area / alpha^(2.0/3.0)'
    pp_names = 'area alpha'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  csv = true
[]
