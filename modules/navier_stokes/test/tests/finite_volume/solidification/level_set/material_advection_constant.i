[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 8
    xmin = 0
    xmax = 1
  []
[]

[Problem]
  linear_sys_names = phi_system
[]

[Variables]
  [level_set_phi]
    type = MooseLinearVariableFVReal
    solver_sys = phi_system
    initial_condition = 3.25
  []
[]

[Functions]
  [expanding_flux]
    type = ParsedFunction
    expression = x
  []
  [exact_constant]
    type = ParsedFunction
    expression = 3.25
  []
[]

[LinearFVKernels]
  [time]
    type = LinearFVTimeDerivative
    variable = level_set_phi
  []
  [material_advection]
    type = LinearFVMaterialAdvection
    variable = level_set_phi
    volumetric_face_flux = expanding_flux
    advected_interp_method = upwind
  []
[]

[LinearFVBCs]
  [outflow]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = level_set_phi
    boundary = right
    use_two_term_expansion = false
  []
[]

[Postprocessors]
  [minimum]
    type = ElementExtremeFunctorValue
    functor = level_set_phi
    value_type = min
    execute_on = FINAL
  []
  [maximum]
    type = ElementExtremeFunctorValue
    functor = level_set_phi
    value_type = max
    execute_on = FINAL
  []
  [l2_error]
    type = ElementL2FunctorError
    approximate = level_set_phi
    exact = exact_constant
    execute_on = FINAL
  []
[]

[Executioner]
  type = Transient
  system_names = phi_system
  scheme = bdf2
  dt = 0.05
  num_steps = 4
  l_tol = 1e-13
[]

[Outputs]
  csv = true
  file_base = material_advection_constant
  execute_on = FINAL
[]
