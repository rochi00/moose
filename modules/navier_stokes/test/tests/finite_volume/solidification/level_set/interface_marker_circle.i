# A stationary signed-distance circle exercises interface-band AMR independently of advection and
# redistancing. The base 8x8 mesh is refined three times only near phi = 0.

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 8
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
[]

[Problem]
  solve = false
  linear_sys_names = transport_system
[]

[Variables]
  [level_set_phi]
    type = MooseLinearVariableFVReal
    solver_sys = transport_system
  []
[]

[Functions]
  [circle]
    type = ParsedFunction
    expression = '0.2 - sqrt((x - 0.5)^2 + (y - 0.5)^2)'
  []
[]

[ICs]
  [circle]
    type = FunctionIC
    variable = level_set_phi
    function = circle
  []
[]

[AuxVariables]
  [h_level]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [h_level]
    type = ElementAdaptivityLevelAux
    variable = h_level
    level = h
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Adaptivity]
  initial_marker = interface_band
  initial_steps = 3
  max_h_level = 3
  [Markers]
    [interface_band]
      type = LevelSetInterfaceMarker
      level_set_variable = level_set_phi
      refine_band_cells = 2
      coarsen_band_cells = 4
    []
  []
[]

[Postprocessors]
  [active_elements]
    type = NumElements
    execute_on = 'INITIAL FINAL'
  []
  [maximum_h_level]
    type = ElementExtremeValue
    variable = h_level
    value_type = max
    execute_on = 'INITIAL FINAL'
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
  exodus = true
  file_base = interface_marker_circle
  execute_on = 'INITIAL FINAL'
[]
