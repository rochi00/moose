# Forward-reverse vortex benchmark from CASL-p4est
# examples/ml_mass_conservation/vortex_2d.cpp.
#
# CASL uses a circle centered at (0.5, 0.75) with radius 0.15, a finest spacing
# of 1/64, CFL = 1, ten reinitialization iterations, and a velocity reversal at
# t = 0.625. The exact interface at t = 1.25 is therefore the initial circle.
# The base 2 x 2 mesh with five refinement levels reproduces CASL's level 1 to
# level 6 mesh hierarchy. Our phi sign is positive inside to match the melting
# model; CASL's sign is negative inside.

base_cells = 2
max_level = 5
radius = 0.15
exact_area = 0.07068583470577035
center_x = 0.5
center_y = 0.75
half_time = 0.625
end_time = 1.25
finest_h = ${fparse 1.0 / base_cells / 2^max_level}
interface_width = ${fparse 1.5 * finest_h}
dt = ${finest_h}

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${base_cells}
    ny = ${base_cells}
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
[]

[Problem]
  linear_sys_names = transport_system
[]

[Variables]
  [level_set_phi]
    type = MooseLinearVariableFVReal
    solver_sys = transport_system
  []
[]

[AuxVariables]
  [material_fraction]
    type = MooseVariableFVReal
  []
  [h_level]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [initial_circle]
    type = ParsedFunction
    expression = '${radius} - sqrt((x - ${center_x})^2 + (y - ${center_y})^2)'
  []
  [vortex_u]
    type = ParsedFunction
    expression = 's:=if(t <= ${half_time}, 1, -1); '
                 '-s * sin(pi*x)^2 * sin(2*pi*y)'
  []
  [vortex_v]
    type = ParsedFunction
    expression = 's:=if(t <= ${half_time}, 1, -1); '
                 's * sin(pi*y)^2 * sin(2*pi*x)'
  []
[]

[ICs]
  [circle]
    type = FunctionIC
    variable = level_set_phi
    function = initial_circle
  []
[]

[FunctorMaterials]
  [properties]
    type = GenericFunctorMaterial
    prop_names = 'material_density zero_flux'
    prop_values = '1 0'
  []
  [regularized_heaviside]
    type = ParsedFunctorMaterial
    property_name = alpha
    functor_names = level_set_phi
    functor_symbols = phi
    expression = 'if(phi <= -${interface_width}, 0, '
                 'if(phi >= ${interface_width}, 1, '
                 '0.5 * (1 + phi / ${interface_width} + '
                 'sin(pi * phi / ${interface_width}) / pi)))'
  []
  [exact_regularized_heaviside]
    type = ParsedFunctorMaterial
    property_name = exact_alpha
    expression = 'd:=${radius} - sqrt((x - ${center_x})^2 + (y - ${center_y})^2); '
                 'if(d <= -${interface_width}, 0, '
                 'if(d >= ${interface_width}, 1, '
                 '0.5 * (1 + d / ${interface_width} + '
                 'sin(pi * d / ${interface_width}) / pi)))'
  []
  [sharp_indicator]
    type = ParsedFunctorMaterial
    property_name = sharp_alpha
    functor_names = level_set_phi
    functor_symbols = phi
    expression = 'if(phi >= 0, 1, 0)'
  []
  [shape_difference]
    type = ParsedFunctorMaterial
    property_name = shape_difference
    functor_names = 'alpha exact_alpha'
    expression = 'abs(alpha - exact_alpha)'
  []
  [moment_x]
    type = ParsedFunctorMaterial
    property_name = moment_x
    functor_names = alpha
    expression = 'x * alpha'
  []
  [moment_y]
    type = ParsedFunctorMaterial
    property_name = moment_y
    functor_names = alpha
    expression = 'y * alpha'
  []
  [masked_level_set]
    type = ParsedFunctorMaterial
    property_name = masked_level_set
    functor_names = level_set_phi
    functor_symbols = phi
    expression = 'd:=${radius} - sqrt((x - ${center_x})^2 + (y - ${center_y})^2); '
                 'if(abs(d) <= 3 * ${finest_h}, phi, d)'
  []
  [band_level_set_absolute_error]
    type = ParsedFunctorMaterial
    property_name = band_level_set_absolute_error
    functor_names = masked_level_set
    functor_symbols = phi
    expression = 'd:=${radius} - sqrt((x - ${center_x})^2 + (y - ${center_y})^2); '
                 'abs(phi - d)'
  []
[]

[LinearFVKernels]
  [time]
    type = LinearFVTimeDerivative
    variable = level_set_phi
  []
  [advection]
    type = LinearFVMaterialAdvection
    variable = level_set_phi
    u = vortex_u
    v = vortex_v
    advected_interp_method = venkatakrishnan
  []
[]

# The vortex has zero normal velocity on every domain boundary, so no advective
# inflow boundary condition is required.

[AuxKernels]
  [material_fraction]
    type = FunctorAux
    variable = material_fraction
    functor = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [h_level]
    type = ElementAdaptivityLevelAux
    variable = h_level
    level = h
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[UserObjects]
  [level_set_correction]
    type = MassConservedLevelSetCorrector
    level_set_variable = level_set_phi
    material_fraction = alpha
    material_density = material_density
    volumetric_face_flux = zero_flux
    interface_width = ${interface_width}
    relative_tolerance = 1e-11
    redistance_iterations = 10
    report = true
    execute_on = TIMESTEP_END
    execution_order_group = -1
  []
[]

[Adaptivity]
  initial_marker = interface_band
  initial_steps = ${max_level}
  marker = interface_band
  cycles_per_step = 1
  max_h_level = ${max_level}
  [Markers]
    [interface_band]
      type = LevelSetInterfaceMarker
      level_set_variable = level_set_phi
      # Two CASL finest-cell diagonals are 2*sqrt(2)*h in two dimensions.
      refine_band_cells = 3
      coarsen_band_cells = 5
      conservative_prolongation = true
      conservative_interface_width = ${interface_width}
      level_set_corrector = level_set_correction
    []
  []
[]

[Postprocessors]
  [active_elements]
    type = NumElements
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [maximum_h_level]
    type = ElementExtremeValue
    variable = h_level
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [material_mass]
    type = ElementIntegralFunctorPostprocessor
    functor = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mass_relative_error]
    type = ChangeOverTimePostprocessor
    postprocessor = material_mass
    change_with_respect_to_initial = true
    compute_relative_change = true
    take_absolute_value = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [maximum_mass_relative_error]
    type = TimeExtremeValue
    postprocessor = mass_relative_error
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sharp_area]
    type = ElementIntegralFunctorPostprocessor
    functor = sharp_alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sharp_area_relative_error]
    type = ParsedPostprocessor
    pp_names = sharp_area
    expression = 'abs(sharp_area - ${exact_area}) / ${exact_area}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sharp_area_relative_change]
    type = ChangeOverTimePostprocessor
    postprocessor = sharp_area
    change_with_respect_to_initial = true
    compute_relative_change = true
    take_absolute_value = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shape_l1_error]
    type = ElementIntegralFunctorPostprocessor
    functor = shape_difference
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [band_level_set_l2_error]
    type = ElementL2FunctorError
    approximate = masked_level_set
    exact = initial_circle
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [band_level_set_linf_error]
    type = ElementExtremeFunctorValue
    functor = band_level_set_absolute_error
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [first_moment_x]
    type = ElementIntegralFunctorPostprocessor
    functor = moment_x
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [first_moment_y]
    type = ElementIntegralFunctorPostprocessor
    functor = moment_y
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [centroid_error]
    type = ParsedPostprocessor
    pp_names = 'first_moment_x first_moment_y material_mass'
    expression = 'sqrt((first_moment_x/material_mass - ${center_x})^2 + '
                 '(first_moment_y/material_mass - ${center_y})^2)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  system_names = transport_system
  dt = ${dt}
  end_time = ${end_time}
  l_tol = 1e-12
  [TimeIntegrator]
    type = BDF2
    # Do not let the first reverse-flow step use pre-reversal BDF2 history.
    restart_times = '${half_time}'
  []
[]

[Outputs]
  csv = true
  exodus = true
  file_base = distorted_circle_casl_p4est
  execute_on = 'INITIAL TIMESTEP_END'
  # This also makes the history restart exact if dt is overridden adaptively.
  sync_times = '${half_time}'
  [console]
    type = Console
    execute_postprocessors_on = 'INITIAL FINAL'
  []
[]
