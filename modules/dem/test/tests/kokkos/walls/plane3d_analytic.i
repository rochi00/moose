# The analytic-plane reference of plane3d.i (V12): the same drop on the plane y = 0 instead of the
# triangulated 'bottom' sideset. See plane3d.i.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 2
    nz = 4
    xmax = 0.4
    ymax = 0.2
    zmax = 0.4
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.08 0.03 0.09'
    initial_velocity = '0.5 -0.3 0.3'
    radius = 0.01
    density = 2500
    gravity = '0 -9.81 0'
    normal_stiffness = 1e5
    normal_damping = 1
    tangential_stiffness = 1e5
    tangential_damping = 1
    friction = 0.3
    wall_points = '0 0 0'
    wall_normals = '0 1 0'
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
  dt = 0.05
  num_steps = 4
[]

[Outputs]
  csv = true
[]
