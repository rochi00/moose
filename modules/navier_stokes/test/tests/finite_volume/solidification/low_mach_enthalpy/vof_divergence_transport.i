[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 16
  []
[]

[Problem]
  linear_sys_names = 'u_system pressure_system alpha_system mass_system energy_system'
  previous_nl_solution_required = true
[]

[Variables]
  [rho_adv]
    type = MooseLinearVariableFVReal
    solver_sys = mass_system
  []
  [T]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
  []
[]

[AuxVariables]
  [divergence_source_aux]
    type = MooseVariableFVReal
  []
  [alpha_divergence_aux]
    type = MooseVariableFVReal
  []
[]

[ICs]
  [rho_adv]
    type = FunctionIC
    variable = rho_adv
    function = initial_density
  []
  [temperature]
    type = ConstantIC
    variable = T
    value = 5
  []
[]

[Functions]
  [material_indicator_initial]
    type = ParsedFunction
    expression = 'if(x < 0.5, 1, 0)'
  []
  [initial_density]
    type = ParsedFunction
    expression = 'if(x < 0.5, 8, 0.1)'
  []
  [localized_heat_source]
    type = ParsedFunction
    # Keep the phase-change divergence away from the gas/PCM interface so bulk preservation and
    # interface transport are independently observable.
    expression = 'if(x < 0.25, 584, 0)'
  []
  [zero_face_flux]
    type = ConstantFunction
    value = 0
  []
  [pressure_initial]
    type = ConstantFunction
    value = 0
  []
[]

[FunctorMaterials]
  [constants]
    type = GenericFunctorMaterial
    prop_names = 'cp_gas cp_solid cp_liquid rho_gas rho_solid rho_liquid k_gas k_solid k_liquid mu_gas mu_solid mu_liquid latent_heat T_solidus T_liquidus T_reference'
    prop_values = '2 3 5 0.1 8 6 0.1 1 0.5 0.01 1 1 100 10 14 0'
  []
  [thermodynamics]
    type = LowMachEnthalpyThermodynamicsFunctorMaterial
    material_fraction = alpha
    temperature = T
    cp_gas = cp_gas
    cp_solid = cp_solid
    cp_liquid = cp_liquid
    rho_gas = rho_gas
    rho_solid = rho_solid
    rho_liquid = rho_liquid
    k_gas = k_gas
    k_solid = k_solid
    k_liquid = k_liquid
    mu_gas = mu_gas
    mu_solid = mu_solid
    mu_liquid = mu_liquid
    L = latent_heat
    T_solidus = T_solidus
    T_liquidus = T_liquidus
    T_reference = T_reference
  []
  [alpha_divergence]
    type = ParsedFunctorMaterial
    # This is the cellwise source in d(alpha)/dt + div(alpha*u) = alpha*div(u).
    property_name = alpha_divergence
    functor_names = 'alpha low_mach_divergence_source'
    expression = 'alpha * low_mach_divergence_source'
  []
[]

[AuxKernels]
  [divergence_source_aux]
    type = FunctorAux
    variable = divergence_source_aux
    functor = low_mach_divergence_source
    execute_on = 'TIMESTEP_END FINAL'
  []
  [alpha_divergence_aux]
    type = FunctorAux
    variable = alpha_divergence_aux
    functor = alpha_divergence
    execute_on = 'TIMESTEP_END FINAL'
  []
[]

[Physics]
  [NavierStokes]
    [ConservativeSharpInterfaceFlowSegregated]
      [flow]
        velocity_variable = vel_x
        pressure_variable = pressure
        compressibility = incompressible
        density = rho_adv
        dynamic_viscosity = low_mach_dynamic_viscosity
        volume_fraction_functor = alpha
        gravity = '0 0 0'
        initial_velocity = '0 0 0'
        initial_pressure = pressure_initial
        wall_boundaries = left
        momentum_wall_types = noslip
        outlet_boundaries = right
        momentum_outlet_types = fixed-pressure
        pressure_functors = 0
        orthogonality_correction = false
        momentum_two_term_bc_expansion = false
        pressure_two_term_bc_expansion = false
        momentum_advection_interpolation = upwind
        momentum_mass_flux = low_mach_mass_flux
      []
    []
    [ConservativeSharpInterfaceVOFSegregated]
      [material_indicator]
        coupled_flow_physics = flow
        passive_scalar_names = alpha
        system_names = alpha_system
        initial_scalar_variables = material_indicator_initial
        liquid_density_name = rho_liquid
        gas_density_name = rho_gas
        liquid_dynamic_viscosity_name = mu_liquid
        gas_dynamic_viscosity_name = mu_gas
        mixture_density_name = vof_mixture_density
        mixture_dynamic_viscosity_name = vof_mixture_dynamic_viscosity
        compression_factor = 1
        interface_normal_functor = flow_interface_unit_normal_face
        phase_change_divergence = low_mach_divergence_source
        n_alpha_corrections = 1
        n_limiter_iterations = 3
        fail_on_unbounded_alpha = true
      []
    []
  []
[]

[LinearFVKernels]
  [mass_time]
    type = LinearFVTimeDerivative
    variable = rho_adv
  []
  [mass_advection]
    type = LinearFVLowMachMassAdvection
    variable = rho_adv
    volumetric_face_flux = vof_transport_phi
  []
  [enthalpy_time]
    type = LinearFVLowMachEnthalpyTimeDerivative
    variable = T
    temporary_density = rho_adv
    density = low_mach_density
    specific_enthalpy = low_mach_enthalpy_from_temperature
    dh_dT = low_mach_dh_dT
  []
  [enthalpy_advection]
    type = LinearFVLowMachEnthalpyAdvection
    variable = T
    specific_enthalpy = low_mach_enthalpy_from_temperature
    dh_dT = low_mach_dh_dT
    mass_flux_functor = low_mach_mass_flux
    inflow_specific_enthalpy = 10
  []
  [heat_source]
    type = LinearFVSource
    variable = T
    source_density = localized_heat_source
  []
  [low_mach_divergence]
    type = LinearFVLowMachDivergence
    variable = pressure
    divergence_source = low_mach_divergence_source
  []
  [solid_drag_x]
    type = LinearFVLowMachSolidDrag
    variable = vel_x
    drag_coefficient = low_mach_solid_drag
  []
[]

[LinearFVBCs]
  [mass_boundary]
    type = LinearFVInletOutletScalarBC
    variable = rho_adv
    boundary = 'left right'
    face_flux = vof_transport_phi
    backflow_value = 0.1
  []
  [temperature_boundary]
    type = LinearFVInletOutletScalarBC
    variable = T
    boundary = 'left right'
    face_flux = zero_face_flux
    backflow_value = 5
  []
[]

[UserObjects]
  [mass_flux]
    type = LowMachImplicitMassFlux
    system = mass_system
    mass_advection_kernel = mass_advection
    temporary_density = rho_adv
    eos_density = low_mach_density
    mass_flux_name = low_mach_mass_flux
    minimum_density = 0.09
    maximum_density = 8
  []
  [enthalpy_state]
    type = LowMachEnthalpyNewtonState
    temperature_variable = T
    material_fraction = alpha
    specific_enthalpy = low_mach_enthalpy_from_temperature
    dh_dT = low_mach_dh_dT
    liquid_fraction = low_mach_liquid_fraction
    cp_gas = cp_gas
    cp_solid = cp_solid
    cp_liquid = cp_liquid
    rho_solid = rho_solid
    rho_liquid = rho_liquid
    h_solid = low_mach_solid_enthalpy
    h_liquid = low_mach_liquid_enthalpy
    T_solidus = T_solidus
    T_liquidus = T_liquidus
    T_reference = T_reference
  []
  [divergence_source]
    type = LowMachDivergenceSource
    system = energy_system
    enthalpy_advection_kernel = enthalpy_advection
    temporary_density = rho_adv
    eos_density = low_mach_density
    specific_enthalpy = low_mach_enthalpy_from_temperature
    dliquid_fraction_dh = low_mach_dliquid_fraction_dh
    material_fraction = alpha
    rho_solid = rho_solid
    rho_liquid = rho_liquid
    h_solid = low_mach_solid_enthalpy
    h_liquid = low_mach_liquid_enthalpy
  []
[]

[Executioner]
  type = ReducedPressurePIMPLE
  scheme = bdf2
  rhie_chow_user_object = ins_rhie_chow_interpolator

  momentum_systems = u_system
  pressure_system = pressure_system
  energy_system = energy_system
  active_scalar_systems = alpha_system
  low_mach_mass_system = mass_system
  low_mach_mass_flux = mass_flux
  low_mach_enthalpy_state = enthalpy_state
  low_mach_divergence_source = divergence_source

  num_iterations = 2
  num_piso_iterations = 0
  continue_on_max_its = true
  startup_pressure_initialization = none
  active_scalar_equation_relaxation = 1

  momentum_absolute_tolerance = 1e-12
  pressure_absolute_tolerance = 1e-12
  energy_absolute_tolerance = 1e-12
  low_mach_mass_absolute_tolerance = 1e-12
  active_scalar_absolute_tolerance = 1e-12

  momentum_l_tol = 1e-12
  pressure_l_tol = 1e-12
  energy_l_tol = 1e-12
  low_mach_mass_l_tol = 1e-12
  active_scalar_l_tol = 1e-12

  momentum_petsc_options_iname = '-pc_type'
  momentum_petsc_options_value = lu
  pressure_petsc_options_iname = '-pc_type'
  pressure_petsc_options_value = lu
  energy_petsc_options_iname = '-pc_type'
  energy_petsc_options_value = lu
  low_mach_mass_petsc_options_iname = '-pc_type'
  low_mach_mass_petsc_options_value = lu
  active_scalar_petsc_options_iname = '-pc_type'
  active_scalar_petsc_options_value = lu

  dt = 1
  num_steps = 1
[]

[Postprocessors]
  [material_fraction]
    type = ElementIntegralVariablePostprocessor
    variable = alpha
  []
  [minimum_material_fraction]
    type = ElementExtremeValue
    variable = alpha
    value_type = min
  []
  [maximum_material_fraction]
    type = ElementExtremeValue
    variable = alpha
    value_type = max
  []
  [left_bulk_material_fraction]
    type = PointValue
    variable = alpha
    point = '0.03125 0 0'
  []
  [right_bulk_material_fraction]
    type = PointValue
    variable = alpha
    point = '0.96875 0 0'
  []
  [temperature]
    type = ElementAverageValue
    variable = T
  []
  [temporary_density]
    type = ElementAverageValue
    variable = rho_adv
  []
  [enthalpy]
    type = ElementIntegralFunctorPostprocessor
    functor = low_mach_enthalpy_from_temperature
  []
  [liquid_fraction]
    type = ElementIntegralFunctorPostprocessor
    functor = low_mach_liquid_fraction
  []
  [eos_density]
    type = ElementIntegralFunctorPostprocessor
    functor = low_mach_density
  []
  [divergence_source_value]
    type = ElementExtremeFunctorValue
    functor = low_mach_divergence_source
  []
  [left_alpha_flux]
    type = SideIntegralFunctorPostprocessor
    boundary = left
    functor = alpha_phi_limited
  []
  [right_alpha_flux]
    type = SideIntegralFunctorPostprocessor
    boundary = right
    functor = alpha_phi_limited
  []
  [integrated_alpha_divergence]
    type = ElementIntegralVariablePostprocessor
    variable = alpha_divergence_aux
  []
  [solid_drag]
    type = ElementExtremeFunctorValue
    functor = low_mach_solid_drag
  []
[]

[Outputs]
  csv = true
[]
