Nx = 12
max_level = 2
finest_h = ${fparse 1.0 / Nx / 2^max_level}
interface_width = 1e-6

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
  solve = false
[]

[Variables]
  [level_set_phi]
    type = MooseLinearVariableFVReal
    solver_sys = transport_system
  []
[]

[AuxVariables]
  [h_level]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [distorted_phi]
    type = ParsedFunction
    expression = 'd:=sqrt((x - 0.5)^2 + (y - 0.5)^2) - 0.25; '
                 '(1 + 0.25 * sin(2*pi*x) * sin(2*pi*y)) * d'
  []
  [exact_phi]
    type = ParsedFunction
    expression = 'sqrt((x - 0.5)^2 + (y - 0.5)^2) - 0.25'
  []
[]

[ICs]
  [level_set]
    type = FunctionIC
    variable = level_set_phi
    function = distorted_phi
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
  [masked_level_set]
    type = ParsedFunctorMaterial
    property_name = masked_level_set
    functor_names = level_set_phi
    functor_symbols = phi
    expression = 'd:=sqrt((x - 0.5)^2 + (y - 0.5)^2) - 0.25; '
                 'if(abs(d) <= 2 * ${finest_h}, phi, d)'
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
    relative_tolerance = 1e-12
    redistance_iterations = 50
    execute_on = TIMESTEP_BEGIN
  []
[]

[Adaptivity]
  initial_marker = interface_band
  initial_steps = ${max_level}
  max_h_level = ${max_level}
  [Markers]
    [interface_band]
      type = LevelSetInterfaceMarker
      level_set_variable = level_set_phi
      refine_band_cells = 3
      coarsen_band_cells = 6
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
  [band_distance_l2_error]
    type = ElementL2FunctorError
    approximate = masked_level_set
    exact = exact_phi
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [material_mass]
    type = ElementIntegralFunctorPostprocessor
    functor = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mass_relative_change]
    type = ChangeOverTimePostprocessor
    postprocessor = material_mass
    change_with_respect_to_initial = true
    compute_relative_change = true
    take_absolute_value = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  system_names = transport_system
  scheme = implicit-euler
  dt = 0.01
  num_steps = 1
[]

[Outputs]
  csv = true
  exodus = true
  file_base = redistance_adaptive_circle
  execute_on = 'INITIAL TIMESTEP_END'
[]
