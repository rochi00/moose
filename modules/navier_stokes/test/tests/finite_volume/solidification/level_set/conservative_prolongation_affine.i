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
    expression = 'x + 2 * y'
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
    []
  []
[]

[Postprocessors]
  [active_elements]
    type = NumElements
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [l2_error]
    type = ElementL2FunctorError
    approximate = level_set_phi
    exact = affine
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [field_integral]
    type = ElementIntegralVariablePostprocessor
    variable = level_set_phi
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
  file_base = conservative_prolongation_affine
  execute_on = 'INITIAL TIMESTEP_END'
[]
