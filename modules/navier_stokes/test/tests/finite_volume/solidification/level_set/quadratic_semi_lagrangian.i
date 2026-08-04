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
  linear_sys_names = transport_system
[]

[Variables]
  [level_set_phi]
    type = MooseLinearVariableFVReal
    solver_sys = transport_system
  []
[]

[Functions]
  [initial_quadratic]
    type = ParsedFunction
    expression = '1 + 2*x - 3*y + 4*x^2 + 5*x*y - 2*y^2'
  []
  [velocity_x]
    type = ParsedFunction
    expression = '0.5'
  []
  [velocity_y]
    type = ParsedFunction
    expression = '-0.25'
  []
  [exact_quadratic]
    type = ParsedFunction
    expression = 'xp:=x-0.5*t; yp:=y+0.25*t; '
                 '1 + 2*xp - 3*yp + 4*xp^2 + 5*xp*yp - 2*yp^2'
  []
[]

[ICs]
  [level_set_phi]
    type = FunctionIC
    variable = level_set_phi
    function = initial_quadratic
  []
[]

[FunctorMaterials]
  [error]
    type = ParsedFunctorMaterial
    property_name = absolute_error
    functor_names = 'level_set_phi exact_quadratic'
    functor_symbols = 'phi exact'
    expression = 'abs(phi-exact)'
  []
[]

[LinearFVKernels]
  [semi_lagrangian]
    type = LinearFVQuadraticSemiLagrangianAdvection
    variable = level_set_phi
    u = velocity_x
    v = velocity_y
  []
[]

[Postprocessors]
  [maximum_error]
    type = ElementExtremeFunctorValue
    functor = absolute_error
    value_type = max
    execute_on = FINAL
  []
[]

[Executioner]
  type = Transient
  system_names = transport_system
  dt = 0.05
  num_steps = 1
  l_tol = 1e-13
[]

[Outputs]
  csv = true
  file_base = quadratic_semi_lagrangian
  execute_on = FINAL
[]
