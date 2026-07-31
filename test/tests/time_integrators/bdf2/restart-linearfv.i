# Linear FV counterpart of restart.i. The source changes discontinuously at
# t = 0.5, so the first step in the second interval must restart BDF2.
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Problem]
  linear_sys_names = u_system
[]

[Variables]
  [u]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
  []
[]

[Functions]
  [piecewise_source]
    type = ParsedFunction
    expression = 'if(t <= 0.5, 1, -1)'
  []
[]

[LinearFVKernels]
  [time]
    type = LinearFVTimeDerivative
    variable = u
  []
  [source]
    type = LinearFVSource
    variable = u
    source_density = piecewise_source
  []
[]

[Postprocessors]
  [numerical_u]
    type = ElementAverageValue
    variable = u
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  system_names = u_system
  dt = 0.25
  end_time = 1
  l_tol = 1e-13
  [TimeIntegrator]
    type = BDF2
    restart_times = '0.5'
  []
[]

[Outputs]
  csv = true
  sync_times = '0.5'
  execute_on = 'INITIAL TIMESTEP_END'
[]
