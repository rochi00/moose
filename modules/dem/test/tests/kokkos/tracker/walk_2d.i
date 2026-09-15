# Particles advected across a 10x10 QUAD4 mesh, crossing 2.5 cells per step in x and 1.5 in y,
# so the face walk must take several hops per step.
# Every step the device face-walk assignment is checked against libMesh's PointLocator
# (verify = true errors on any mismatch).
# Particles crossing a partition boundary are migrated to the owning rank, so the output does
# not depend on the number of ranks. Particles reaching the right or top edge are counted as lost.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [tracker]
    type = KokkosParticleTracker
    initial_positions = '0.05 0.05 0
                         0.13 0.47 0
                         0.31 0.29 0
                         0.52 0.71 0
                         0.05 0.95 0
                         0.66 0.12 0
                         0.85 0.55 0
                         0.97 0.03 0'
    velocity = '2.5 1.5 0'
    verify = true
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [num_particles]
    type = KokkosParticleTrackerValue
    tracker = tracker
    value = num_particles
    execute_on = TIMESTEP_END
  []
  [num_lost]
    type = KokkosParticleTrackerValue
    tracker = tracker
    value = num_lost
    execute_on = TIMESTEP_END
  []
  [max_hops]
    type = KokkosParticleTrackerValue
    tracker = tracker
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
