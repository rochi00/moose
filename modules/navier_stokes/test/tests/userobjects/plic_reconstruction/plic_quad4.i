[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx = 2
    ny = 1
  []
[]

[Problem]
  linear_sys_names = 'alpha_system phi_system'
  previous_nl_solution_required = true
  solve = false
[]

[Variables]
  [alpha]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_system
  []
  [phi]
    type = MooseLinearVariableFVReal
    solver_sys = phi_system
  []
[]

[Functions]
  [alpha_init]
    type = ParsedFunction
    expression = 'if(x < 0.5, 1, 0.4)'
  []
  [phi_init]
    type = ParsedFunction
    expression = 'x'
  []
[]

[ICs]
  [alpha_ic]
    type = FunctionIC
    variable = alpha
    function = alpha_init
  []
  [phi_ic]
    type = FunctionIC
    variable = phi
    function = phi_init
  []
[]

[UserObjects]
  [plic]
    type = PLICReconstruction
    alpha_variable = alpha
    level_set_variable = phi
    execute_on = initial
  []
[]

[Postprocessors]
  [plane_error]
    type = PLICReconstructionVolumeError
    plic_uo = plic
    alpha_variable = alpha
    execute_on = initial
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
[]
