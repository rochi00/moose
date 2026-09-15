# V1: a freely spinning sphere conserves its spin angular momentum and rotational kinetic energy.
# With a scalar moment of inertia and no torque the angular velocity is constant, so the gold is
# the analytic  I = 2/5 m r^2,  KE_rot = 1/2 I |w|^2,  |L| = I |w|.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.5 0.5 0'
    initial_angular_velocity = '1 2 3'
    radius = 0.1
    density = 1000
    substeps = 10000
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [rotational_kinetic_energy]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = rotational_kinetic_energy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [angular_momentum]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = angular_momentum
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  dt = 1
  num_steps = 2
[]

[Outputs]
  csv = true
[]
