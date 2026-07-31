interface_width = 0.025

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 80
    xmin = 0
    xmax = 1
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
  [velocity]
    type = ConstantFunction
    value = 1
  []
  [exact_phi]
    type = ParsedFunction
    expression = '0.1 - abs(x - 0.35 - t)'
  []
[]

[ICs]
  [level_set]
    type = FunctionIC
    variable = level_set_phi
    function = exact_phi
  []
[]

[FunctorMaterials]
  [properties]
    type = GenericFunctorMaterial
    prop_names = material_density
    prop_values = 1
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

[LinearFVKernels]
  [time]
    type = LinearFVTimeDerivative
    variable = level_set_phi
  []
  [advection]
    type = LinearFVMaterialAdvection
    variable = level_set_phi
    u = velocity
    advected_interp_method = venkatakrishnan
  []
[]

[LinearFVBCs]
  [inflow]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = level_set_phi
    boundary = left
    functor = exact_phi
  []
  [outflow]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = level_set_phi
    boundary = right
    use_two_term_expansion = false
  []
[]

[UserObjects]
  [mass_correction]
    type = MassConservedLevelSetCorrector
    level_set_variable = level_set_phi
    material_fraction = alpha
    material_density = material_density
    volumetric_face_flux = velocity
    interface_width = ${interface_width}
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [material_mass]
    type = ElementIntegralFunctorPostprocessor
    functor = alpha
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [mass_relative_change]
    type = ChangeOverTimePostprocessor
    postprocessor = material_mass
    change_with_respect_to_initial = true
    compute_relative_change = true
    take_absolute_value = true
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [alpha_minimum]
    type = ElementExtremeFunctorValue
    functor = alpha
    value_type = min
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [alpha_maximum]
    type = ElementExtremeFunctorValue
    functor = alpha
    value_type = max
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Executioner]
  type = Transient
  system_names = transport_system
  scheme = implicit-euler
  dt = 0.005
  num_steps = 2
  l_tol = 1e-13
[]

[Outputs]
  csv = true
  file_base = mass_conserved_translation
  execute_on = 'INITIAL TIMESTEP_END'
[]
