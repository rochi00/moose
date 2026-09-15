# A sphere launched sliding along a wall with Coulomb friction: the friction force mu m g slows it
# and its torque spins it up until the contact point stops slipping at t* = 2 v0 / (7 mu g) = 0.097 s,
# after which it rolls at 5/7 v0 (generate.py). The sphere starts at rest on the linear spring, at
# the static overlap m g / k, so the normal contact is quiet throughout; the tangential spring is
# saturated at mu m g while slipping and holds the contact point at rest afterwards, the tangential
# dashpot damping out the stick oscillation.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 1
    xmax = 0.5
    ymax = 0.1
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.1 0.0099989725 0'
    initial_velocity = '1 0 0'
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
    substeps = 1500
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
  dt = 0.03
  num_steps = 10
[]

[Outputs]
  csv = true
[]
