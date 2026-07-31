interface_location = 0.35
epsilon = 0.05

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 40
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
  [constant_velocity]
    type = ParsedFunction
    expression = 1
  []
  [exact_phi]
    type = ParsedFunction
    expression = 'x - ${interface_location} - t'
  []
[]

[ICs]
  [phi]
    type = FunctionIC
    variable = level_set_phi
    function = exact_phi
  []
[]

[FunctorMaterials]
  [regularized_heaviside]
    type = ParsedFunctorMaterial
    property_name = pcm_heaviside
    functor_names = level_set_phi
    functor_symbols = phi
    expression = 'if(phi < -${epsilon}, 0, if(phi > ${epsilon}, 1, '
                 '0.5 * (1 + phi / ${epsilon} + sin(pi * phi / ${epsilon}) / pi)))'
  []
[]

[LinearFVKernels]
  [phi_time]
    type = LinearFVTimeDerivative
    variable = level_set_phi
  []
  [phi_advection]
    type = LinearFVMaterialAdvection
    variable = level_set_phi
    u = constant_velocity
    advected_interp_method = upwind
  []
[]

[LinearFVBCs]
  [phi_inflow]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = level_set_phi
    boundary = left
    functor = exact_phi
  []
  [phi_outflow]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = level_set_phi
    boundary = right
    use_two_term_expansion = false
  []
[]

[Postprocessors]
  [phi_l2_error]
    type = ElementL2FunctorError
    approximate = level_set_phi
    exact = exact_phi
    execute_on = FINAL
  []
  [heaviside_minimum]
    type = ElementExtremeFunctorValue
    functor = pcm_heaviside
    value_type = min
    execute_on = FINAL
  []
  [heaviside_maximum]
    type = ElementExtremeFunctorValue
    functor = pcm_heaviside
    value_type = max
    execute_on = FINAL
  []
[]

[Executioner]
  type = Transient
  system_names = transport_system
  scheme = bdf2
  dt = 0.025
  num_steps = 4
  l_tol = 1e-13
[]

[Outputs]
  csv = true
  file_base = material_advection_translation
  execute_on = FINAL
[]
