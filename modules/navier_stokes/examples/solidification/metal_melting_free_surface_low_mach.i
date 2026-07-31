# Section 6.2 metal-melting benchmark from:
# R. Thirumalaisamy and A. P. S. Bhalla, Int. J. Multiphase Flow 169 (2023) 104605.
#
# The paper uses a reinitialized level set for the gas/metal surface. This branch currently uses
# its conservative bounded VOF transport for the same gas/metal material indicator. The
# unequal-density enthalpy closure, phase-change divergence, temporary-density mass equation, and
# common mass flux follow the low-Mach implementation.
#
# MOOSE finite-volume variables do not currently accept periodic constraints. The paper's data are
# uniform in x, so zero-normal-flux symmetry sides give the same one-dimensional solution. They
# must be replaced by true periodic sides if lateral perturbations are added in a future FV path.
#
# Paper-resolution defaults require a long run: 256 x 256 cells, dt = 1e-3 s, t_end = 250 s.
# A quick wiring check can be run with, for example:
#   navier_stokes-opt -i metal_melting_free_surface_low_mach.i \
#     Nx=2 Ny=20 Executioner/num_steps=1 Outputs/exodus/sync_only=false

Nx = 256
Ny = 256
dt = 1e-3
end_time = 250

rho_gas = 0.4
rho_solid = 2475
rho_liquid = 2700

cp_gas = 1100
cp_solid = 910
cp_liquid = 1042.4

k_gas = 6.1e-2
k_solid = 211
k_liquid = 91

mu_gas = 4e-5
mu_solid = 1.4e-3
mu_liquid = 1.4e-3

latent_heat = 383840
T_reference = 933.6
T_solidus = 928.6
T_liquidus = 938.6
T_cold = 840.24
T_hot_liquid = 4668
T_bottom = 5601.6

initial_liquid_height = 0.30
initial_metal_height = 0.45
expected_final_height = 0.4375
cell_height = ${fparse 1.0 / Ny}
gas_inflow_specific_enthalpy = ${fparse cp_gas * (T_cold - T_reference)}

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx = ${Nx}
    ny = ${Ny}
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system alpha_system mass_system energy_system'
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
  [eos_density]
    type = MooseVariableFVReal
  []
  [phase_change_divergence]
    type = MooseVariableFVReal
  []
  [pcm_mass_density]
    type = MooseVariableFVReal
  []
  [pcm_liquid_volume_density]
    type = MooseVariableFVReal
  []
  [pcm_solid_volume_density]
    type = MooseVariableFVReal
  []
[]

[Functions]
  [initial_material_fraction]
    type = ParsedFunction
    # Exact cell-volume fraction for the initially horizontal free surface.
    expression = 'min(max((metal_height - (y - 0.5 * dy)) / dy, 0), 1)'
    symbol_names = 'metal_height dy'
    symbol_values = '${initial_metal_height} ${cell_height}'
  []
  [initial_temperature]
    type = ParsedFunction
    expression = 'if(y <= liquid_height, T_liquid, T_other)'
    symbol_names = 'liquid_height T_liquid T_other'
    symbol_values = '${initial_liquid_height} ${T_hot_liquid} ${T_cold}'
  []
  [initial_density]
    type = ParsedFunction
    expression = 'rho_gas + (if(y <= liquid_height, rho_liquid, rho_solid) - rho_gas) * min(max((metal_height - (y - 0.5 * dy)) / dy, 0), 1)'
    symbol_names = 'liquid_height metal_height dy rho_liquid rho_solid rho_gas'
    symbol_values = '${initial_liquid_height} ${initial_metal_height} ${cell_height} ${rho_liquid} ${rho_solid} ${rho_gas}'
  []
  [pressure_initial]
    type = ConstantFunction
    value = 0
  []
[]

[ICs]
  [rho_adv]
    type = FunctionIC
    variable = rho_adv
    function = initial_density
  []
  [temperature]
    type = FunctionIC
    variable = T
    function = initial_temperature
  []
[]

[FunctorMaterials]
  [constants]
    type = GenericFunctorMaterial
    prop_names = 'cp_gas cp_solid cp_liquid rho_gas rho_solid rho_liquid k_gas k_solid k_liquid mu_gas mu_solid mu_liquid latent_heat T_solidus T_liquidus T_reference'
    prop_values = '${cp_gas} ${cp_solid} ${cp_liquid} ${rho_gas} ${rho_solid} ${rho_liquid} ${k_gas} ${k_solid} ${k_liquid} ${mu_gas} ${mu_solid} ${mu_liquid} ${latent_heat} ${T_solidus} ${T_liquidus} ${T_reference}'
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
  [pcm_mass_density]
    type = ParsedFunctorMaterial
    property_name = pcm_mass_density_functor
    functor_names = 'alpha low_mach_liquid_fraction'
    functor_symbols = 'H fl'
    expression = 'H * (${rho_solid} + (${rho_liquid} - ${rho_solid}) * fl)'
  []
  [pcm_phase_density]
    type = ParsedFunctorMaterial
    property_name = pcm_phase_density
    functor_names = low_mach_liquid_fraction
    functor_symbols = fl
    expression = '${rho_solid} + (${rho_liquid} - ${rho_solid}) * fl'
  []
  [pcm_liquid_volume_density]
    type = ParsedFunctorMaterial
    property_name = pcm_liquid_volume_density_functor
    functor_names = 'alpha low_mach_liquid_fraction'
    functor_symbols = 'H fl'
    expression = 'H * fl'
  []
  [pcm_solid_volume_density]
    type = ParsedFunctorMaterial
    property_name = pcm_solid_volume_density_functor
    functor_names = 'alpha low_mach_liquid_fraction'
    functor_symbols = 'H fl'
    expression = 'H * (1 - fl)'
  []
[]

[Physics]
  [NavierStokes]
    [ConservativeSharpInterfaceFlowSegregated]
      [flow]
        velocity_variable = 'vel_x vel_y'
        pressure_variable = pressure
        compressibility = incompressible
        density = rho_adv
        dynamic_viscosity = low_mach_dynamic_viscosity
        volume_fraction_functor = alpha
        gravity = '0 0 0'
        initial_velocity = '0 0 0'
        initial_pressure = pressure_initial

        wall_boundaries = 'left right bottom'
        momentum_wall_types = 'symmetry symmetry noslip'
        outlet_boundaries = top
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
        initial_scalar_variables = initial_material_fraction
        volume_fraction_outlet_type = inlet-outlet

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
        n_limiter_iterations = 6
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
    material_compression_flux = alpha_compression_flux
    material_density = pcm_phase_density
    background_density = rho_gas
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
    inflow_specific_enthalpy = ${gas_inflow_specific_enthalpy}
  []
  [thermal_diffusion]
    type = LinearFVDiffusion
    variable = T
    diffusion_coeff = low_mach_thermal_conductivity
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
  [solid_drag_y]
    type = LinearFVLowMachSolidDrag
    variable = vel_y
    drag_coefficient = low_mach_solid_drag
  []
[]

[LinearFVBCs]
  [mass_top]
    type = LinearFVInletOutletScalarBC
    variable = rho_adv
    boundary = top
    face_flux = vof_transport_phi
    backflow_value = ${rho_gas}
  []
  [temperature_bottom]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = T
    boundary = bottom
    functor = ${T_bottom}
  []
  [temperature_top]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = T
    boundary = top
  []
[]

[AuxKernels]
  [liquid_fraction]
    type = FunctorAux
    variable = liquid_fraction
    functor = low_mach_liquid_fraction
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [eos_density]
    type = FunctorAux
    variable = eos_density
    functor = low_mach_density
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [phase_change_divergence]
    type = FunctorAux
    variable = phase_change_divergence
    functor = low_mach_divergence_source
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pcm_mass_density]
    type = FunctorAux
    variable = pcm_mass_density
    functor = pcm_mass_density_functor
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pcm_liquid_volume_density]
    type = FunctorAux
    variable = pcm_liquid_volume_density
    functor = pcm_liquid_volume_density_functor
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pcm_solid_volume_density]
    type = FunctorAux
    variable = pcm_solid_volume_density
    functor = pcm_solid_volume_density_functor
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
    minimum_density = 0.39
    maximum_density = ${rho_liquid}
    density_tolerance = 1e-3
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

  momentum_systems = 'u_system v_system'
  pressure_system = pressure_system
  energy_system = energy_system
  active_scalar_systems = alpha_system
  low_mach_mass_system = mass_system
  low_mach_mass_flux = mass_flux
  low_mach_enthalpy_state = enthalpy_state
  low_mach_divergence_source = divergence_source

  low_mach_enthalpy_max_iterations = 20
  low_mach_enthalpy_tolerance = 1e-8
  low_mach_enthalpy_residual_tolerance = 1e-8

  num_iterations = 12
  num_piso_iterations = 0
  continue_on_max_its = false
  startup_pressure_initialization = none
  active_scalar_equation_relaxation = 1

  momentum_absolute_tolerance = 1e-8
  pressure_absolute_tolerance = 1e-8
  energy_absolute_tolerance = 1e-8
  low_mach_mass_absolute_tolerance = 1e-8
  active_scalar_absolute_tolerance = 1e-8

  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  energy_l_abs_tol = 1e-12
  low_mach_mass_l_abs_tol = 1e-12
  active_scalar_l_abs_tol = 1e-12

  momentum_l_tol = 1e-10
  pressure_l_tol = 1e-10
  energy_l_tol = 1e-10
  low_mach_mass_l_tol = 1e-10
  active_scalar_l_tol = 1e-10

  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
  low_mach_mass_petsc_options_iname = '-pc_type -pc_hypre_type'
  low_mach_mass_petsc_options_value = 'hypre boomeramg'
  active_scalar_petsc_options_iname = '-pc_type -pc_hypre_type'
  active_scalar_petsc_options_value = 'hypre boomeramg'

  dt = ${dt}
  end_time = ${end_time}
[]

[Postprocessors]
  [pcm_volume]
    type = ElementIntegralVariablePostprocessor
    variable = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pcm_mass]
    type = ElementIntegralVariablePostprocessor
    variable = pcm_mass_density
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pcm_liquid_volume]
    type = ElementIntegralVariablePostprocessor
    variable = pcm_liquid_volume_density
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pcm_solid_volume]
    type = ElementIntegralVariablePostprocessor
    variable = pcm_solid_volume_density
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pcm_mass_error_percent]
    type = ParsedPostprocessor
    pp_names = pcm_mass_relative_change
    expression = 'pcm_mass_relative_change * 100'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pcm_mass_relative_change]
    type = ChangeOverTimePostprocessor
    postprocessor = pcm_mass
    change_with_respect_to_initial = true
    compute_relative_change = true
    take_absolute_value = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [final_height_error]
    type = ParsedPostprocessor
    pp_names = pcm_volume
    expression = 'pcm_volume - ${expected_final_height}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [minimum_material_fraction]
    type = ElementExtremeValue
    variable = alpha
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [maximum_material_fraction]
    type = ElementExtremeValue
    variable = alpha
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [maximum_velocity]
    type = ElementExtremeValue
    variable = vel_y
    value_type = max_abs
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [maximum_divergence_source]
    type = ElementExtremeValue
    variable = phase_change_divergence
    value_type = max_abs
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [minimum_temporary_density]
    type = ElementExtremeValue
    variable = rho_adv
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [maximum_temporary_density]
    type = ElementExtremeValue
    variable = rho_adv
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [exodus]
    type = Exodus
    execute_on = 'INITIAL TIMESTEP_END'
    sync_times = '100 150 250'
    sync_only = true
  []
[]
