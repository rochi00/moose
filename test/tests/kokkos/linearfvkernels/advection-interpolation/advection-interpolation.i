# Advection-dominated advection-diffusion of a smoothed step on a skewed mesh with each advected
# interpolation method, selected with the top-level variable interp. Each time step reassembles
# the steady system with the weights and deferred corrections of the previous solution (there is
# no time derivative kernel). The gold files come from this input with the host objects: the types
# of the kernels, the boundary conditions and the functions replaced by LinearFVAdvection,
# LinearFVDiffusion, LinearFVAdvectionDiffusionFunctorDirichletBC,
# LinearFVAdvectionDiffusionFunctorNeumannBC and ParsedFunction.

interp = upwind

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 12
    ny = 10
  []
  # Smooth perturbation of the interior nodes that vanishes on the boundary
  [skewed]
    type = ParsedNodeTransformGenerator
    input = gmg
    x_function = 'x + 0.03 * sin(3.141592653589793 * x) * sin(6.283185307179586 * y)'
    y_function = 'y + 0.03 * sin(6.283185307179586 * x) * sin(3.141592653589793 * y)'
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
  []
[]

[ICs]
  [u_ic]
    type = FunctionIC
    variable = u
    function = 'x * y'
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
    limit_to_linear = false
    blending_factor = 0.8
  []
  [muscl]
    type = FVAdvectedMUSCLDeferredCorrection
    deferred_correction_factor = 0.9
  []
[]

[LinearFVKernels]
  [advection]
    type = KokkosLinearFVAdvection
    variable = u
    velocity = '1 0.6 0'
    advected_interp_method_name = ${interp}
  []
  [diffusion]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = diffusivity
  []
[]

[LinearFVBCs]
  [inlet]
    type = KokkosLinearFVFunctorDirichletBC
    variable = u
    boundary = 'left bottom'
    functor = inlet_value
  []
  [outlet]
    type = KokkosLinearFVFunctorNeumannBC
    variable = u
    boundary = 'right top'
    functor = zero
    diffusion_coeff = diffusivity
  []
[]

[Functions]
  [diffusivity]
    type = KokkosParsedFunction
    expression = '0.01'
  []
  [zero]
    type = KokkosParsedFunction
    expression = '0'
  []
  [inlet_value]
    type = KokkosParsedFunction
    expression = '1 / (1 + exp(-40 * (y - 0.4)))'
  []
[]

[Executioner]
  type = Transient
  system_names = u_sys
  num_steps = 4
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  file_base = advection-interpolation_${interp}
  exodus = true
[]
