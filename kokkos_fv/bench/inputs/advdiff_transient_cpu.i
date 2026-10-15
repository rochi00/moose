# Transient CPU reference for advdiff_transient_kokkos.i (CPU linear FV objects, same settings as
# advdiff_cpu.i plus LinearFVTimeDerivative with implicit Euler).

n = 16
elem_type = HEX8
dt = 0.01
skew = 0.0
steps = 50
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
    type = LinearFVTimeDerivative
    variable = u
  []
  [advection]
    type = LinearFVAdvection
    variable = u
    velocity = '1 0.5 0.25'
    advected_interp_method_name = ${interp}
  []
  [diffusion]
    type = LinearFVDiffusion
    variable = u
    diffusion_coeff = k_fn
    use_nonorthogonal_correction = ${noc}
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
