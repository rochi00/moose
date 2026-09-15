# V8: a sphere rolling without slipping on a horizontal plane at 1 m/s with rolling resistance.
# The rolling spring saturates at once (the contact point is at rest, so there is no tangential
# transient) and the plastic torque mu_r m g r decelerates it at mu_r g / 1.4 = 0.70 m/s^2 until
# it stops after 0.714 m at t = 1.43 s, where the spring holds it at rest (generate.py). The
# sphere starts at the static overlap m g / k on the linear spring.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 1
    xmax = 1
    ymax = 0.1
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.1 0.009998972699202276 0'
    initial_velocity = '1 0 0'
    initial_angular_velocity = '0 0 -100'
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
    rolling_friction = 0.1
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
  dt = 0.1
  num_steps = 20
[]

[Outputs]
  csv = true
[]
