# Properties of Kokkos linear FV objects given as numbers, KokkosConstantFunction and
# KokkosParsedFunction objects, with the diffusion coefficient interpolated to the faces with the
# harmonic and the geometric average, on a skewed mesh with two block-restricted variables. The
# boundary conditions on the interface between the blocks use functions that are discontinuous
# there, which tells where they are evaluated: as on the host, the Neumann flux of u from the side
# of u (1), the Dirichlet value of v given by a function without a side, at the interface (2), and
# given by a functor material without a side, as the interpolation of its cell values (between 1
# and 2). The gold file comes from this input with the host
# objects: the types of the kernels, the boundary conditions and the functions replaced by
# LinearFVDiffusion, LinearFVSource, LinearFVAdvectionDiffusionFunctorDirichletBC,
# LinearFVAdvectionDiffusionFunctorNeumannBC, ParsedFunctorMaterial, ParsedFunction and
# ConstantFunction. The tests run it with each interface value of v.

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '0.5 0.5'
    ix = '6 6'
    dy = '1'
    iy = '10'
    subdomain_id = '1 2'
  []
  # Smooth perturbation of the nodes that vanishes on the boundary and on the interface
  [skewed]
    type = ParsedNodeTransformGenerator
    input = cmg
    x_function = 'x + 0.02 * sin(6.283185307179586 * x) * sin(6.283185307179586 * y)'
    y_function = 'y + 0.02 * sin(6.283185307179586 * x) * sin(3.141592653589793 * y)'
    z_function = 'z'
  []
  # One sideset on the interface, whose boundary conditions apply to the variables on both sides
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = skewed
    new_boundary = interface
    primary_block = 1
    paired_block = 2
  []
[]

[Problem]
  linear_sys_names = 'u_sys'
[]

[Variables]
  [u]
    type = MooseLinearVariableFVReal
    solver_sys = 'u_sys'
    block = 1
  []
  [v]
    type = MooseLinearVariableFVReal
    solver_sys = 'u_sys'
    block = 2
  []
[]

[FVInterpolationMethods]
  [average]
    type = FVGeometricAverage
  []
  [harmonic]
    type = FVHarmonicAverage
  []
[]

[FunctorMaterials]
  [step]
    type = KokkosParsedFunctorMaterial
    property_name = step_material
    expression = 'if(x > 0.5, 1, 2)'
  []
[]

[LinearFVKernels]
  [diffusion_u]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = k_fn
    coeff_interp_method = harmonic
    block = 1
  []
  [source_u]
    type = KokkosLinearFVSource
    variable = u
    source_density = s_fn
    scaling_factor = 2
    block = 1
  []
  [diffusion_v]
    type = KokkosLinearFVDiffusion
    variable = v
    diffusion_coeff = k_fn
    coeff_interp_method = average
    block = 2
  []
  [source_v]
    type = KokkosLinearFVSource
    variable = v
    source_density = 3
    scaling_factor = half
    block = 2
  []
[]

[LinearFVBCs]
  [outer_u]
    type = KokkosLinearFVFunctorDirichletBC
    variable = u
    boundary = left
    functor = 0
  []
  [outer_v]
    type = KokkosLinearFVFunctorDirichletBC
    variable = v
    boundary = right
    functor = 1
  []
  [interface_u]
    type = KokkosLinearFVFunctorNeumannBC
    variable = u
    boundary = interface
    functor = step_flux
    diffusion_coeff = k_fn
  []
  [interface_v]
    type = KokkosLinearFVFunctorDirichletBC
    variable = v
    boundary = interface
    functor = step_value
  []
  [flux_u]
    type = KokkosLinearFVFunctorNeumannBC
    variable = u
    boundary = 'top bottom'
    functor = flux
    diffusion_coeff = k_fn
  []
  [flux_v]
    type = KokkosLinearFVFunctorNeumannBC
    variable = v
    boundary = 'top bottom'
    functor = flux
    diffusion_coeff = k_fn
  []
[]

[Functions]
  [k_fn]
    type = KokkosParsedFunction
    expression = '1 + 4 * x + y * y'
  []
  [s_fn]
    type = KokkosParsedFunction
    expression = '1 + x * y'
  []
  [half]
    type = KokkosConstantFunction
    value = 0.5
  []
  [step_flux]
    type = KokkosParsedFunction
    expression = 'if(x < 0.5, 1, -1)'
  []
  [step_value]
    type = KokkosParsedFunction
    expression = 'if(x > 0.5, 1, 2)'
  []
  [flux]
    type = KokkosParsedFunction
    expression = '0.1 * x'
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
