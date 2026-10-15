# Green-Gauss gradient of a linear FV variable computed with Kokkos, on a non-uniform skewed mesh
# with Dirichlet boundary conditions on two sides, a Neumann boundary condition on one and none on
# the last. The diffusion kernel uses the non-orthogonal correction, which reads the gradients. The
# gold file comes from this input with the host objects: the types of the gradient method, the
# kernels, the boundary conditions and the functions replaced by FVGreenGaussGradient,
# LinearFVDiffusion, LinearFVSource, LinearFVAdvectionDiffusionFunctorDirichletBC,
# LinearFVAdvectionDiffusionFunctorNeumannBC and ParsedFunction.
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 6
    bias_x = 1.2
    bias_y = 0.85
  []
  # Smooth perturbation of the interior nodes that vanishes on the boundary, so boundaries stay
  # straight while the interior faces become non-orthogonal
  [skewed]
    type = ParsedNodeTransformGenerator
    input = gmg
    x_function = 'x + 0.04 * sin(3.141592653589793 * x) * sin(6.283185307179586 * y)'
    y_function = 'y + 0.04 * sin(6.283185307179586 * x) * sin(3.141592653589793 * y)'
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
    gradient_method = gg
  []
[]

[FVGradientMethods]
  [gg]
    type = KokkosFVGreenGaussGradient
  []
[]

[AuxVariables]
  [grad_u_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [grad_u_y]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [grad_u_x]
    type = LinearFVGradientComponentTestAux
    variable = grad_u_x
    linear_variable = u
    component = 0
    execute_on = timestep_end
  []
  [grad_u_y]
    type = LinearFVGradientComponentTestAux
    variable = grad_u_y
    linear_variable = u
    component = 1
    execute_on = timestep_end
  []
[]

[LinearFVKernels]
  [diffusion]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = one
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
  [neumann]
    type = KokkosLinearFVFunctorNeumannBC
    variable = u
    boundary = 'top'
    functor = top_flux
    diffusion_coeff = one
  []
[]

[Functions]
  [one]
    type = KokkosParsedFunction
    expression = '1'
  []
  [source]
    type = KokkosParsedFunction
    expression = '4 * x * y'
  []
  [boundary_value]
    type = KokkosParsedFunction
    expression = '1 + 2 * x + y'
  []
  [top_flux]
    type = KokkosParsedFunction
    expression = '1 + x'
  []
[]

[Executioner]
  type = Steady
  system_names = u_sys
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  exodus = true
[]
