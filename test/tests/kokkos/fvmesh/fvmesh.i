# Checks the Kokkos finite volume mesh data against the MOOSE face and element info on a skewed
# mesh with two blocks and a different sideset on each side of their interface
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 3
    nz = 3
  []
  [block]
    type = SubdomainBoundingBoxGenerator
    input = gmg
    bottom_left = '0.5 0 0'
    top_right = '1 1 1'
    block_id = 1
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = block
    primary_block = 0
    paired_block = 1
    new_boundary = interface
  []
  [interface_other_side]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface
    primary_block = 1
    paired_block = 0
    new_boundary = interface_other_side
  []
  [skew]
    type = ParsedNodeTransformGenerator
    input = interface_other_side
    # Affine, so that the faces stay planar and the faces of a cell close its surface exactly
    x_function = 'x + 0.2 * y + 0.1 * z'
    y_function = 'y + 0.1 * x'
    z_function = 'z + 0.15 * x - 0.1 * y'
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

[LinearFVKernels]
  [diffusion]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = one
  []
[]

[LinearFVBCs]
  [dirichlet]
    type = KokkosLinearFVFunctorDirichletBC
    variable = u
    boundary = 'left right'
    functor = one
  []
[]

[Functions]
  [one]
    type = KokkosParsedFunction
    expression = 1
  []
[]

[Postprocessors]
  [failures]
    type = KokkosFVMeshTest
  []
[]

[Executioner]
  type = Steady
  system_names = u_sys
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
