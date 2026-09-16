# Particles are inserted at 400 per second at random non-overlapping positions in a box near the
# top of a 2D domain, fall under gravity with friction onto the floor sideset, and the ones that
# reach the outflow box at the right end of the floor are removed. The counts of particles
# inserted (num_particles + num_exited) follow the rate exactly, no two particles overlap at
# insertion, and the inserted ones settle into a bed against the walls.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
    xmax = 0.2
    ymax = 0.2
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.1 0.01 0'
    radius = 0.005
    density = 2500
    gravity = '0 -9.81 0'
    normal_stiffness = 1e5
    normal_damping = 1
    tangential_stiffness = 1e5
    tangential_damping = 1
    friction = 0.3
    wall_boundaries = 'bottom left right'
    insertion_box = '0.02 0.15 0
                     0.12 0.19 0'
    insertion_rate = 400
    insertion_velocity = '0.5 0 0'
    insertion_seed = 3
    outflow_box = '0.17 0 0
                   0.2 0.02 0'
    skin = 0.001
    substeps = 400
    verify = true
    verify_pair_states = true
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [num_particles]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = num_particles
  []
  [num_exited]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = num_exited
  []
  [num_contacts]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = num_contacts
  []
  [floor]
    type = KokkosWallForce
    cloud = cloud
    boundary = bottom
    component = y
  []
[]

[Executioner]
  type = Transient
  dt = 0.01
  num_steps = 60
[]

[Outputs]
  csv = true
[]
