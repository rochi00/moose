# Forward-reverse deforming-vortex benchmark from the CASL p4est level-set examples.
# The divergence-free flow deforms the disk until half_time, then reverses exactly.
# At end_time the exact interface is the initial circle.

Nx = 64
radius = 0.15
exact_area = 0.07068583470577035
center_x = 0.5
center_y = 0.75
half_time = 0.625
end_time = 1.25
h = ${fparse 1.0 / Nx}
interface_width = ${fparse 1.5 * h}
dt = ${h}

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${Nx}
    ny = ${Nx}
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
  [density]
    type = GenericFunctorMaterial
    prop_names = material_density
    prop_values = 1
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
                 'if(abs(d) <= 3 * ${h}, phi, d)'
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

[UserObjects]
  [level_set_correction]
    type = MassConservedLevelSetCorrector
    level_set_variable = level_set_phi
    material_fraction = alpha
    material_density = material_density
    # The standalone transport kernel constructs flux from vortex_u and vortex_v. This name is
    # metadata for coupling checks used by the segregated flow executioner and is not evaluated by
    # the corrector.
    volumetric_face_flux = vortex_u
    interface_width = ${interface_width}
    relative_tolerance = 1e-11
    redistance_iterations = 10
    report = true
    execute_on = TIMESTEP_END
    execution_order_group = -1
  []
[]

[Postprocessors]
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
  scheme = bdf2
  dt = ${dt}
  end_time = ${end_time}
  l_tol = 1e-12
[]

[Outputs]
  csv = true
  exodus = true
  file_base = distorted_disk_vortex
  execute_on = 'INITIAL TIMESTEP_END'
  [console]
    type = Console
    execute_postprocessors_on = 'INITIAL FINAL'
  []
[]
