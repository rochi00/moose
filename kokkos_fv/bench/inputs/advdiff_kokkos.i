# 3D steady advection-diffusion-source benchmark using the Kokkos linear FV objects.
# Every Transient step reassembles and resolves the same steady system (there is no time
# derivative kernel), which gives a repeated-assembly loop for timing.
# Top-level variables are overridden from the command line by bench.py.

n = 16             # cells per direction
elem_type = HEX8   # HEX8 | TET4 | PRISM6
skew = 0.0         # amplitude of a smooth interior node perturbation (0 = orthogonal mesh)
steps = 10         # number of assemblies/solves

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${n}
    ny = ${n}
    nz = ${n}
    elem_type = ${elem_type}
  []
  # Smooth perturbation that vanishes on the boundary, so boundaries stay planar
  [skewed]
    type = ParsedNodeTransformGenerator
    input = gmg
    x_function = 'x + ${skew} * sin(3.141592653589793 * x) * sin(6.283185307179586 * y) * sin(3.141592653589793 * z) / ${n}'
    y_function = 'y + ${skew} * sin(6.283185307179586 * x) * sin(3.141592653589793 * y) * sin(3.141592653589793 * z) / ${n}'
    z_function = 'z'
  []
[]

[Problem]
  linear_sys_names = 'u_sys'
[]

[Variables]
  [u]
    type = MooseLinearVariableFVReal
    solver_sys = 'u_sys'
    initial_condition = 1.0
  []
[]

[LinearFVKernels]
  [advection]
    type = KokkosLinearFVAdvection
    variable = u
    velocity = '1 0.5 0.25'
  []
  [diffusion]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = k_fn
  []
  [source]
    type = KokkosLinearFVSource
    variable = u
    source_density = s_fn
    # Must be set explicitly: the default "1" builds a KokkosConstantFunction, which the kernel
    # rejects because it requires a KokkosParsedFunction
    scaling_factor = unit_fn
  []
[]

[LinearFVBCs]
  [inflow]
    type = KokkosLinearFVFunctorDirichletBC
    variable = u
    boundary = 'left bottom back'
    functor = g_fn
  []
  [outflow]
    type = KokkosLinearFVFunctorNeumannBC
    variable = u
    boundary = 'right top front'
    functor = zero_fn
  []
[]

[Functions]
  [k_fn]
    type = KokkosParsedFunction
    expression = '1 + 0.5 * x'
  []
  [s_fn]
    type = KokkosParsedFunction
    expression = '1 + x * y * z'
  []
  [g_fn]
    type = KokkosParsedFunction
    expression = '1 + x + y + z'
  []
  [unit_fn]
    type = KokkosParsedFunction
    expression = '1'
  []
  [zero_fn]
    type = KokkosParsedFunction
    expression = '0'
  []
[]

[Executioner]
  type = Transient
  num_steps = ${steps}
  dt = 1
  system_names = u_sys
  # Assembly timing: one unmonitored Richardson/Jacobi sweep keeps the solve negligible and
  # always "converged" (KSP_CONVERGED_ITS when no norm is computed)
  l_max_its = 1
  petsc_options_iname = '-ksp_type -ksp_norm_type -pc_type'
  petsc_options_value = 'richardson none jacobi'
[]

[Postprocessors]
  [ncells]
    type = NumElements
    execute_on = FINAL
  []
[]

[Reporters]
  [perf]
    type = PerfGraphReporter
    execute_on = FINAL
  []
[]

[Outputs]
  console = false
  [json]
    type = JSON
    execute_on = 'FINAL'
  []
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
