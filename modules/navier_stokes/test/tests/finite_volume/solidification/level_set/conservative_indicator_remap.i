interface_width = 0.12

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 4
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
  [plane]
    type = ParsedFunction
    expression = 'x + 0.2 * y - 0.48'
  []
[]

[ICs]
  [level_set]
    type = FunctionIC
    variable = level_set_phi
    function = plane
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
    expression = 'if(phi <= -${interface_width}, 0, if(phi >= ${interface_width}, 1, '
                 '0.5 * (1 + phi / ${interface_width} + '
                 'sin(pi * phi / ${interface_width}) / pi)))'
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
    redistance_iterations = 10
    execute_on = NONE
  []
[]

[LinearFVKernels]
  [time]
    type = LinearFVTimeDerivative
    variable = level_set_phi
  []
[]

[Adaptivity]
  marker = interface_band
  cycles_per_step = 1
  max_h_level = 1
  stop_time = 0.015
  [Markers]
    [interface_band]
      type = LevelSetInterfaceMarker
      level_set_variable = level_set_phi
      refine_band_cells = 100
      coarsen_band_cells = 200
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
  [material_fraction_integral]
    type = ElementIntegralFunctorPostprocessor
    functor = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [material_fraction_minimum]
    type = ElementExtremeFunctorValue
    functor = alpha
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [material_fraction_maximum]
    type = ElementExtremeFunctorValue
    functor = alpha
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [material_fraction_error]
    type = ParsedPostprocessor
    expression = 'abs(material_fraction_integral - 0.62)'
    pp_names = material_fraction_integral
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [level_set_l2_error]
    type = ElementL2FunctorError
    approximate = level_set_phi
    exact = plane
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  system_names = transport_system
  scheme = bdf2
  dt = 0.01
  num_steps = 2
[]

[Outputs]
  csv = true
  file_base = conservative_indicator_remap
  execute_on = 'INITIAL TIMESTEP_END'
[]
