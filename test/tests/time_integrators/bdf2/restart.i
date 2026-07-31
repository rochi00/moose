# The forcing changes discontinuously at t = 0.5. BDF2 must not use the
# pre-discontinuity state when it starts integrating the second time interval.
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Variables]
  [u]
    family = SCALAR
    order = FIRST
    initial_condition = 0
  []
[]

[ScalarKernels]
  [time]
    type = ODETimeDerivative
    variable = u
  []
  [piecewise_source]
    type = ParsedODEKernel
    variable = u
    expression = '-if(source_time <= 0.5, 1, -1)'
    postprocessors = source_time
  []
[]

[Postprocessors]
  [source_time]
    type = TimePostprocessor
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [numerical_u]
    type = ScalarVariable
    variable = u
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  dt = 0.25
  end_time = 1
  solve_type = NEWTON
  nl_abs_tol = 1e-14
  [TimeIntegrator]
    type = BDF2
    restart_times = '0.5'
  []
[]

[Outputs]
  csv = true
  hide = source_time
  sync_times = '0.5'
  execute_on = 'INITIAL TIMESTEP_END'
[]
