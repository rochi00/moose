# Properties of Kokkos linear FV objects given by linear FV variables: the diffusion coefficient of
# u is a KokkosParsedFunctorMaterial of the variable T of another system (with Dirichlet and
# Neumann boundary conditions, whose boundary values enter the coefficient on the boundary faces),
# of an auxiliary linear FV variable, of a number and of space and time, and the source density of
# u is T itself, on a skewed mesh. The gold file comes from this input with the host objects: the
# types of the kernels, the boundary conditions, the functor material and the Kokkos functions
# replaced by LinearFVDiffusion, LinearFVSource, LinearFVAdvectionDiffusionFunctorDirichletBC,
# LinearFVAdvectionDiffusionFunctorNeumannBC, ParsedFunctorMaterial and ParsedFunction.

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
  linear_sys_names = 'T_sys u_sys'
[]

[Variables]
  [T]
    type = MooseLinearVariableFVReal
    solver_sys = 'T_sys'
  []
  [u]
    type = MooseLinearVariableFVReal
    solver_sys = 'u_sys'
  []
[]

[AuxVariables]
  [mu_t]
    type = MooseLinearVariableFVReal
  []
[]

[ICs]
  [mu_t]
    type = FunctionIC
    variable = mu_t
    function = mu_t_ic
  []
[]

[FunctorMaterials]
  [k]
    type = KokkosParsedFunctorMaterial
    property_name = k_u
    expression = 'k0 * (1 + 0.5 * temperature * temperature) + mu_t + 0.1 * x * (1 + t)'
    functor_names = '2 T mu_t'
    functor_symbols = 'k0 temperature mu_t'
  []
[]

[LinearFVKernels]
  [diffusion_T]
    type = KokkosLinearFVDiffusion
    variable = T
    diffusion_coeff = 1
  []
  [source_T]
    type = KokkosLinearFVSource
    variable = T
    source_density = 4
  []
  [diffusion_u]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = k_u
  []
  [source_u]
    type = KokkosLinearFVSource
    variable = u
    source_density = T
  []
[]

[LinearFVBCs]
  [T_dirichlet]
    type = KokkosLinearFVFunctorDirichletBC
    variable = T
    boundary = 'left right'
    functor = T_boundary
  []
  [T_neumann]
    type = KokkosLinearFVFunctorNeumannBC
    variable = T
    boundary = 'top bottom'
    functor = T_flux
  []
  [u_dirichlet]
    type = KokkosLinearFVFunctorDirichletBC
    variable = u
    boundary = 'left right top bottom'
    functor = 0
  []
[]

[Functions]
  [mu_t_ic]
    type = ParsedFunction
    expression = '0.2 + x * y'
  []
  [T_boundary]
    type = KokkosParsedFunction
    expression = '1 + x'
  []
  [T_flux]
    type = KokkosParsedFunction
    expression = '0.5 * x'
  []
[]

[Executioner]
  type = Steady
  system_names = 'T_sys u_sys'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  exodus = true
[]
