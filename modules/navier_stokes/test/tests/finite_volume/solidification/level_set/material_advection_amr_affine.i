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
  [affine]
    type = ParsedFunction
    expression = '1 + x + 2 * y'
  []
  [velocity_x]
    type = ParsedFunction
    expression = '1'
  []
  [velocity_y]
    type = ParsedFunction
    expression = '0.5'
  []
  [manufactured_source]
    type = ParsedFunction
    expression = '3 + x + 2 * y'
  []
[]

[ICs]
  [level_set]
    type = FunctionIC
    variable = level_set_phi
    function = affine
  []
[]

[LinearFVKernels]
  [advection]
    type = LinearFVMaterialAdvection
    variable = level_set_phi
    u = velocity_x
    v = velocity_y
    advected_interp_method = venkatakrishnan
  []
  [source]
    type = LinearFVSource
    variable = level_set_phi
    source_density = manufactured_source
  []
  [reaction]
    type = LinearFVReaction
    variable = level_set_phi
  []
[]

[LinearFVBCs]
  [exact]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = level_set_phi
    boundary = 'left right bottom top'
    functor = affine
  []
[]

[Adaptivity]
  initial_marker = refine_right
  initial_steps = 1
  max_h_level = 1
  [Markers]
    [refine_right]
      type = BoxMarker
      bottom_left = '0.5 0 0'
      top_right = '1 1 0'
      inside = refine
      outside = do_nothing
    []
  []
[]

[Postprocessors]
  [active_elements]
    type = NumElements
    execute_on = FINAL
  []
  [l2_error]
    type = ElementL2FunctorError
    approximate = level_set_phi
    exact = affine
    execute_on = FINAL
  []
  [field_integral]
    type = ElementIntegralVariablePostprocessor
    variable = level_set_phi
    execute_on = FINAL
  []
[]

[Executioner]
  type = Steady
  system_names = transport_system
  l_tol = 1e-13
[]

[Outputs]
  csv = true
  file_base = material_advection_amr_affine
  execute_on = FINAL
[]
