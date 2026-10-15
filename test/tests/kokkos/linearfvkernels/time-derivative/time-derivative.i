# Transient diffusion with a source, a time derivative with a constant multiplier, and Dirichlet
# boundaries. The gold files come from this input with the host linear FV objects: the types of
# the kernels, the boundary condition and the functions replaced by LinearFVTimeDerivative,
# LinearFVDiffusion, LinearFVSource, LinearFVAdvectionDiffusionFunctorDirichletBC and ParsedFunction.
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 6
  []
[]

[Problem]
  linear_sys_names = 'u_sys'
[]

[Variables]
  [u]
    type = MooseLinearVariableFVReal
    solver_sys = 'u_sys'
    initial_condition = 1
  []
[]

[LinearFVKernels]
  [time]
    type = KokkosLinearFVTimeDerivative
    variable = u
    factor = 2
  []
  [diffusion]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = coeff
  []
  [source]
    type = KokkosLinearFVSource
    variable = u
    source_density = source
    scaling_factor = one
  []
[]

[LinearFVBCs]
  [dirichlet]
    type = KokkosLinearFVFunctorDirichletBC
    variable = u
    boundary = 'left right'
    functor = boundary_value
  []
[]

[Functions]
  [coeff]
    type = KokkosParsedFunction
    expression = '1 + x * y'
  []
  [source]
    type = KokkosParsedFunction
    expression = '10 * x * (1 - x)'
  []
  [one]
    type = KokkosParsedFunction
    expression = '1'
  []
  [boundary_value]
    type = KokkosParsedFunction
    expression = '1 + y'
  []
[]

[Executioner]
  type = Transient
  num_steps = 4
  dt = 0.05
  scheme = implicit-euler
  system_names = u_sys
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  exodus = true
[]
