# CPU (host) reference for advdiff_kokkos.i: same mesh, coefficients, and boundary conditions,
# using the standard linear FV objects. Settings chosen to match the Kokkos discretization:
#   - upwind advection
#   - no non-orthogonal correction (the Kokkos diffusion kernel has none)
#   - zero-flux Neumann outflow (the Kokkos Neumann functor is a normal gradient, the CPU one a
#     diffusive flux divided by k; they coincide only for zero flux)

n = 16
elem_type = HEX8
skew = 0.0
steps = 10

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${n}
    ny = ${n}
    nz = ${n}
    elem_type = ${elem_type}
  []
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

[FVInterpolationMethods]
  [upwind]
    type = FVAdvectedUpwind
  []
[]

[LinearFVKernels]
  [advection]
    type = LinearFVAdvection
    variable = u
    velocity = '1 0.5 0.25'
    advected_interp_method_name = upwind
  []
  [diffusion]
    type = LinearFVDiffusion
    variable = u
    diffusion_coeff = k_fn
    use_nonorthogonal_correction = false
  []
  [source]
    type = LinearFVSource
    variable = u
    source_density = s_fn
  []
[]

[LinearFVBCs]
  [inflow]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = u
    boundary = 'left bottom back'
    functor = g_fn
  []
  [outflow]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = u
    boundary = 'right top front'
    functor = 0
    diffusion_coeff = k_fn
  []
[]

[Functions]
  [k_fn]
    type = ParsedFunction
    expression = '1 + 0.5 * x'
  []
  [s_fn]
    type = ParsedFunction
    expression = '1 + x * y * z'
  []
  [g_fn]
    type = ParsedFunction
    expression = '1 + x + y + z'
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
