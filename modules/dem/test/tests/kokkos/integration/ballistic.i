# V2: a particle under gravity with no contact follows the analytic parabola
#   x(t) = 0.1 + 2 t,  y(t) = 0.9 + t - 5 t^2
# Velocity Verlet is exact for constant acceleration, so the gold is the analytic solution and the
# comparison tolerance is roundoff.

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
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.1 0.9 0'
    initial_velocity = '2 1 0'
    radius = 0.01
    density = 1000
    gravity = '0 -10 0'
    substeps = 40
    verify = true
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
  dt = 0.04
  num_steps = 5
[]

[Outputs]
  csv = true
[]
