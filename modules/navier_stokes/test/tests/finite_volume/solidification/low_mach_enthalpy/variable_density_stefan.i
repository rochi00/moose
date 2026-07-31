# Fast one-dimensional validation configuration. The paper-resolution study is
# obtained with Nx=1280, length=1, end_time=10, and the case-specific time step.
Nx = 40
length = 0.02
end_time = 0.1
time_step = 1e-3
sample_start = ${fparse length / (2 * Nx)}
sample_end = ${fparse length - sample_start}
inflow_specific_enthalpy = ${fparse 910 * (928.6 - 933.6) + 0.5 * (910 + 1042.4) * (938.6 - 928.6) + 383840 + 1042.4 * (973.6 - 938.6)}

rho_solid = 2475
rho_liquid = 2475

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 1
    nx = ${Nx}
    xmin = 0
    xmax = ${length}
  []
[]

[Problem]
  linear_sys_names = 'u_system pressure_system mass_system energy_system'
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
  [liquid_fraction]
    type = MooseVariableFVReal
  []
[]

[FVInterpolationMethods]
  [harmonic]
    type = FVHarmonicAverage
  []
[]

[ICs]
  [rho_adv]
    type = ConstantIC
    variable = rho_adv
    value = ${rho_liquid}
  []
  [temperature]
    type = ConstantIC
    variable = T
    value = 973.6
  []
[]

[Functions]
  [pressure_initial]
    type = ConstantFunction
    value = 0
  []
[]

[FunctorMaterials]
  [constants]
    type = GenericFunctorMaterial
    prop_names = 'alpha cp_gas cp_solid cp_liquid rho_gas rho_solid rho_liquid k_gas k_solid k_liquid mu_gas mu_solid mu_liquid latent_heat T_solidus T_liquidus T_reference inflow_specific_enthalpy'
    prop_values = '1 1042.4 910 1042.4 ${rho_liquid} ${rho_solid} ${rho_liquid} 91 211 91 0 0 0 383840 928.6 938.6 933.6 ${inflow_specific_enthalpy}'
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
        vof_rho_phi_functor = ''
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
    rhie_chow_user_object = ins_rhie_chow_interpolator
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
    inflow_specific_enthalpy = inflow_specific_enthalpy
    mass_flux_functor = low_mach_mass_flux
  []
  [thermal_diffusion]
    type = LinearFVDiffusion
    variable = T
    diffusion_coeff = low_mach_thermal_conductivity
    coeff_interp_method = harmonic
    use_nonorthogonal_correction = false
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
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = rho_adv
    boundary = 'left right'
    functor = ${rho_liquid}
  []
  [cold_wall]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = T
    boundary = left
    functor = 298.6
  []
  [far_field]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = T
    boundary = right
  []
[]

[AuxKernels]
  [liquid_fraction]
    type = FunctorAux
    variable = liquid_fraction
    functor = low_mach_liquid_fraction
    execute_on = 'INITIAL TIMESTEP_END'
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
    minimum_density = 1
    maximum_density = 3000
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
  low_mach_mass_system = mass_system
  low_mach_mass_flux = mass_flux
  low_mach_enthalpy_state = enthalpy_state
  low_mach_divergence_source = divergence_source
  low_mach_enthalpy_max_iterations = 20
  low_mach_enthalpy_l_tol = 1e-12

  num_iterations = 2
  num_piso_iterations = 0
  continue_on_max_its = true
  startup_pressure_initialization = none

  momentum_absolute_tolerance = 1e-10
  pressure_absolute_tolerance = 1e-10
  energy_absolute_tolerance = 1e-10
  low_mach_mass_absolute_tolerance = 1e-10

  momentum_l_tol = 1e-10
  pressure_l_tol = 1e-10
  energy_l_tol = 1e-10
  low_mach_mass_l_tol = 1e-10

  momentum_petsc_options_iname = '-pc_type'
  momentum_petsc_options_value = lu
  pressure_petsc_options_iname = '-pc_type'
  pressure_petsc_options_value = lu
  energy_petsc_options_iname = '-pc_type'
  energy_petsc_options_value = lu
  low_mach_mass_petsc_options_iname = '-pc_type'
  low_mach_mass_petsc_options_value = lu

  dt = ${time_step}
  end_time = ${end_time}
[]

[VectorPostprocessors]
  [temperature_profile]
    type = LineValueSampler
    start_point = '${sample_start} 0 0'
    end_point = '${sample_end} 0 0'
    num_points = ${Nx}
    variable = T
    sort_by = x
    execute_on = FINAL
  []
  [liquid_fraction_profile]
    type = LineValueSampler
    start_point = '${sample_start} 0 0'
    end_point = '${sample_end} 0 0'
    num_points = ${Nx}
    variable = liquid_fraction
    sort_by = x
    execute_on = FINAL
  []
  [velocity_profile]
    type = LineValueSampler
    start_point = '${sample_start} 0 0'
    end_point = '${sample_end} 0 0'
    num_points = ${Nx}
    variable = vel_x
    sort_by = x
    execute_on = FINAL
  []
[]

[Postprocessors]
  [maximum_velocity]
    type = ElementExtremeValue
    variable = vel_x
    value_type = max_abs
    execute_on = FINAL
  []
[]

[Outputs]
  exodus = false
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]
