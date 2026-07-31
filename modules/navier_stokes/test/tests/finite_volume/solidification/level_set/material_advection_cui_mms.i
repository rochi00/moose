nx = 20
h = ${fparse 1.0 / nx}

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 1
    nx = ${nx}
    xmin = 0
    xmax = 1
  []
  [interior]
    type = ParsedSubdomainMeshGenerator
    input = mesh
    combinatorial_geometry = 'x > ${h} & x < 1 - ${h}'
    block_id = 1
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
    type = ParsedFunction
    expression = '1'
  []
  [exact_cell_phi]
    type = ParsedFunction
    expression = '(exp(x + ${h} / 2) - exp(x - ${h} / 2)) / ${h}'
  []
  [exact_face_phi]
    type = ParsedFunction
    expression = 'exp(x)'
  []
  [source]
    type = ParsedFunction
    expression = '2 * (exp(x + ${h} / 2) - exp(x - ${h} / 2)) / ${h}'
  []
[]

[ICs]
  [phi]
    type = FunctionIC
    variable = level_set_phi
    function = exact_cell_phi
  []
[]

[LinearFVKernels]
  [phi_advection]
    type = LinearFVMaterialAdvection
    variable = level_set_phi
    u = velocity
    advected_interp_method = quick
  []
  [phi_reaction]
    type = LinearFVReaction
    variable = level_set_phi
  []
  [phi_source]
    type = LinearFVSource
    variable = level_set_phi
    source_density = source
  []
[]

[LinearFVBCs]
  [phi_exact]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = level_set_phi
    boundary = 'left right'
    functor = exact_face_phi
  []
[]

[Postprocessors]
  [phi_l2_error]
    type = ElementL2FunctorError
    approximate = level_set_phi
    exact = exact_cell_phi
    block = 1
    execute_on = FINAL
  []
[]

[Executioner]
  type = Steady
  system_names = transport_system
  l_tol = 1e-13
  fixed_point_min_its = 2
  fixed_point_max_its = 30
  fixed_point_force_norms = true
  fixed_point_rel_tol = 1e-10
  fixed_point_abs_tol = 1e-11
[]

[Outputs]
  csv = true
  file_base = material_advection_cui_mms
  execute_on = FINAL
[]
