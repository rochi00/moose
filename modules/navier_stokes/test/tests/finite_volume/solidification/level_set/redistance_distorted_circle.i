Nx = 48
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

[Functions]
  [distorted_phi]
    type = ParsedFunction
    expression = '((x - 0.5)^2 + (y - 0.5)^2 + 0.1) * '
                 '(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 0.25)'
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
                 'if(abs(d) <= 1.5 / ${Nx}, phi, d)'
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
    execute_on = TIMESTEP_BEGIN
  []
[]

[Postprocessors]
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
  l_tol = 1e-13
[]

[Outputs]
  csv = true
  file_base = redistance_distorted_circle
  execute_on = 'INITIAL TIMESTEP_END'
[]
