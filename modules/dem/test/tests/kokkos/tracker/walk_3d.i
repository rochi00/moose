# Three-dimensional version of walk_2d: particles cross 2.5, 1.5, and 2 cells per step in x, y,
# and z through a 10x10x10 HEX8 mesh, so the face walk, migration across 3D partitions, and
# removal on exit are exercised in every direction. Verified against the point locator every step.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 10
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.05 0.05 0.05
                         0.13 0.47 0.22
                         0.31 0.29 0.61
                         0.52 0.71 0.08
                         0.05 0.95 0.44
                         0.66 0.12 0.37
                         0.85 0.55 0.19
                         0.97 0.03 0.90
                         0.42 0.42 0.95
                         0.20 0.80 0.50'
    initial_velocity = '2.5 1.5 2'
    radius = 0.01
    density = 1000
    verify = true
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [num_particles]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = num_particles
    execute_on = TIMESTEP_END
  []
  [num_exited]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = num_exited
    execute_on = TIMESTEP_END
  []
  [num_unresolved]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = num_unresolved
    execute_on = TIMESTEP_END
  []
  [max_hops]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = max_hops
    execute_on = TIMESTEP_END
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  num_steps = 4
[]

[Outputs]
  csv = true
[]
