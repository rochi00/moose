# V14: a sphere resting on the 'bottom' sideset of a 2D mesh under gravity, the floor driven by
# the prescribed nodal displacements disp_x = 0.02 t and disp_y = 0.01 t. The sphere starts at
# the static overlap with the floor's velocity, so it simply rides along: its velocity stays
# (0.02, 0.01), its spin zero (rolling resistance locks it; without the moving floor it would
# fly off and fall back under gravity), and its position advances with the floor to roundoff.
# The wall_reach covers the total displacement so the faces stay candidates of the elements.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 4
    xmax = 0.2
    ymax = 0.2
  []
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxKernels]
  [disp_x]
    type = FunctionAux
    variable = disp_x
    function = '0.02 * t'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [disp_y]
    type = FunctionAux
    variable = disp_y
    function = '0.01 * t'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.05 0.009998972699202276 0'
    initial_velocity = '0.02 0.01 0'
    radius = 0.01
    density = 2500
    gravity = '0 -9.81 0'
    normal_stiffness = 1e5
    normal_damping = 1
    tangential_stiffness = 1e5
    tangential_damping = 1
    friction = 0.5
    rolling_stiffness = 1e5
    rolling_damping = 1
    rolling_friction = 1
    wall_boundaries = 'bottom'
    wall_displacements = 'disp_x disp_y'
    wall_reach = 0.02
    skin = 0.001
    substeps = 5000
    verify = true
    verify_pair_states = true
    execute_on = TIMESTEP_END
  []
[]

[VectorPostprocessors]
  [state]
    type = KokkosParticleState
    cloud = cloud
    execute_on = TIMESTEP_END
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  num_steps = 5
[]

[Outputs]
  csv = true
[]
