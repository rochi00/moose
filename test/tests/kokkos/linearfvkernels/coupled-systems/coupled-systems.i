# Two linear systems solved one after the other: u diffuses, and v diffuses with a source equal to
# u, which the Kokkos kernel of v reads from the other system. The gold file comes from this input
# with the host linear FV objects: the types of the kernels, the boundary conditions and the
# functions replaced by LinearFVDiffusion, LinearFVSource (with source_density = u),
# LinearFVAdvectionDiffusionFunctorDirichletBC and ParsedFunction, run with -w for the unused
# 'coupled' parameter.
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 6
  []
[]

[Problem]
  linear_sys_names = 'u_sys v_sys'
[]

[Variables]
  [u]
    type = MooseLinearVariableFVReal
    solver_sys = 'u_sys'
  []
  [v]
    type = MooseLinearVariableFVReal
    solver_sys = 'v_sys'
  []
[]

[LinearFVKernels]
  [u_diffusion]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = one
  []
  [v_diffusion]
    type = KokkosLinearFVDiffusion
    variable = v
    diffusion_coeff = one
  []
  [v_source]
    type = KokkosLinearFVCoupledSourceTest
    variable = v
    coupled = u
  []
[]

[LinearFVBCs]
  [u]
    type = KokkosLinearFVFunctorDirichletBC
    variable = u
    boundary = 'left right'
    functor = u_boundary
  []
  [v]
    type = KokkosLinearFVFunctorDirichletBC
    variable = v
    boundary = 'left right'
    functor = zero
  []
[]

[Functions]
  [one]
    type = KokkosParsedFunction
    expression = '1'
  []
  [zero]
    type = KokkosParsedFunction
    expression = '0'
  []
  [u_boundary]
    type = KokkosParsedFunction
    expression = '1 + 2 * x + y'
  []
[]

[Executioner]
  type = Steady
  system_names = 'u_sys v_sys'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  exodus = true
[]
