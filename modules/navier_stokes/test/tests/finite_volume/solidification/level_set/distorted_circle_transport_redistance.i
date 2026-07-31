!include distorted_circle_transport_only.i

# Fixed uniform mesh with the production redistancing and global mass projection enabled.
# This separates correction error from mesh-transfer error.

[FunctorMaterials]
  [correction_properties]
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
    report = false
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
  [shape_l1_error]
    type = ElementIntegralFunctorPostprocessor
    functor = shape_difference
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
