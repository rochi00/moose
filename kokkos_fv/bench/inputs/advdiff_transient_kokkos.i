# Transient 3D advection-diffusion-source with implicit Euler, Kokkos objects.
#   du/dt + div(v u) - div(k grad u) = s,  u(t=0) = 1
# The time derivative is KokkosLinearFVTimeDerivative. For the "mixed" variant bench.py swaps it
# for the CPU LinearFVTimeDerivative via LinearFVKernels/time/type=LinearFVTimeDerivative.
# Each step is solved to l_tol with GMRES + Jacobi.

n = 16             # cells per direction
elem_type = HEX8   # HEX8 | TET4 | PRISM6
dt = 0.01
skew = 0.0         # amplitude of a smooth interior node perturbation (0 = orthogonal mesh)
steps = 50         # number of time steps
noc = false        # non-orthogonal correction in the diffusion kernel
interp = upwind    # advected interpolation: upwind | average | minmod | vanleer | muscl

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

[FVInterpolationMethods]
  [upwind]
    type = FVAdvectedUpwind
  []
  [average]
    type = FVGeometricAverage
  []
  [minmod]
    type = FVAdvectedMinmodWeightBased
  []
  [vanleer]
    type = FVAdvectedVanLeerWeightBased
  []
  [muscl]
    type = FVAdvectedMUSCLDeferredCorrection
  []
[]

[LinearFVKernels]
  [time]
    type = KokkosLinearFVTimeDerivative
    variable = u
  []
  [advection]
    type = KokkosLinearFVAdvection
    variable = u
    velocity = '1 0.5 0.25'
    advected_interp_method_name = ${interp}
  []
  [diffusion]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = k_fn
    use_nonorthogonal_correction = ${noc}
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
    diffusion_coeff = k_fn
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
  dt = ${dt}
  system_names = u_sys
  l_tol = 1e-10
  l_max_its = 500
  petsc_options_iname = '-ksp_type -pc_type'
  petsc_options_value = 'gmres jacobi'
[]

[Postprocessors]
  [ncells]
    type = NumElements
    execute_on = FINAL
  []
  [u_avg]
    type = ElementAverageFunctorPostprocessor
    functor = u
    execute_on = TIMESTEP_END
  []
  [u_max]
    type = ElementExtremeFunctorValue
    functor = u
    execute_on = TIMESTEP_END
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
    execute_on = 'TIMESTEP_END FINAL'
  []
[]
