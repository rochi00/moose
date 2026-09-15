# The analytic-wall reference of hopper.i: the two planes through the apex of the V, with inward
# normals (1, 1) / sqrt(2) and (-1, 1) / sqrt(2), instead of the rotated square's sidesets.

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 8
    xmin = -0.5
    xmax = 0.5
    ymin = -0.5
    ymax = 0.5
  []
  [rotate]
    type = TransformGenerator
    input = square
    transform = ROTATE
    vector_value = '45 0 0'
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '-0.12  0.05 0
                          0.10  0.12 0
                         -0.03  0.25 0
                          0.16  0.30 0
                         -0.20  0.38 0
                          0.04  0.45 0'
    radius = 0.03
    density = 2500
    gravity = '0 -9.81 0'
    normal_stiffness = 1e5
    normal_damping = 20
    tangential_stiffness = 1e5
    tangential_damping = 20
    friction = 0.4
    rolling_stiffness = 1e5
    rolling_damping = 20
    rolling_friction = 0.1
    wall_points = '0 -0.7071067811865476 0
                   0 -0.7071067811865476 0'
    wall_normals = '1 1 0
                    -1 1 0'
    skin = 0.003
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
  num_steps = 30
[]

[Outputs]
  csv = true
[]
