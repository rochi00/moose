# Advection-diffusion on a skewed mesh with the Kokkos Robin, outflow, scalar symmetry and
# extrapolated boundary conditions. two_term selects the one-term or two-term expansions, and the
# active list the Dirichlet or extrapolated bottom boundary condition. Each time step reassembles the steady system with
# the gradients and cell values of the previous solution, which the explicit terms of the boundary
# conditions read (there is no time derivative kernel). The gold files come from this input with
# the host objects: the types of the kernels, the boundary conditions and the Kokkos functions
# replaced by LinearFVAdvection, LinearFVDiffusion, LinearFVSource,
# LinearFVAdvectionDiffusionFunctorRobinBC, LinearFVAdvectionDiffusionOutflowBC,
# LinearFVAdvectionDiffusionScalarSymmetryBC, LinearFVAdvectionDiffusionExtrapolatedBC,
# LinearFVAdvectionDiffusionFunctorDirichletBC and ParsedFunction.

two_term = false

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 8
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
    function = u_ic
  []
[]

[FVInterpolationMethods]
  [average]
    type = FVGeometricAverage
  []
[]

[LinearFVKernels]
  [advection]
    type = KokkosLinearFVAdvection
    variable = u
    velocity = '1 0.5 0'
    advected_interp_method_name = average
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
  []
[]

[LinearFVBCs]
  [robin]
    type = KokkosLinearFVFunctorRobinBC
    variable = u
    boundary = left
    alpha = k_fn
    beta = 2
    gamma = gamma_fn
  []
  [outflow]
    type = KokkosLinearFVOutflowBC
    variable = u
    boundary = right
    use_two_term_expansion = ${two_term}
  []
  [symmetry]
    type = KokkosLinearFVScalarSymmetryBC
    variable = u
    boundary = top
    use_two_term_expansion = ${two_term}
  []
  [dirichlet]
    type = KokkosLinearFVFunctorDirichletBC
    variable = u
    boundary = bottom
    functor = '1'
  []
  [extrapolated]
    type = KokkosLinearFVExtrapolatedBC
    variable = u
    boundary = bottom
    use_two_term_expansion = ${two_term}
  []
  active = 'robin outflow symmetry dirichlet'
[]

[Functions]
  [k_fn]
    type = KokkosParsedFunction
    expression = '0.2 + 0.1 * x * y'
  []
  [s_fn]
    type = KokkosParsedFunction
    expression = '1 + x'
  []
  [gamma_fn]
    type = KokkosParsedFunction
    expression = '1 + y * y'
  []
  [u_ic]
    type = ParsedFunction
    expression = '1 + x * y'
  []
[]

[Executioner]
  type = Transient
  system_names = u_sys
  num_steps = 3
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  exodus = true
[]
