# A sphere resting on the floor (an analytic plane at y = 0) is compressed by a servo-controlled
# plane above it driven to a target force of 1 N with the linear spring-dashpot model. The servo
# moves the top wall at a velocity proportional to the force error, capped, and stops when the
# measured force reaches the target: at rest the sphere is squeezed between the walls with the
# top overlap F / k and the floor overlap (F + m g) / k, so the top wall ends at
# y = 2 r - (F + m g) / k - F / k from the floor. Gravity adds the weight to the floor's force.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
    xmax = 0.1
    ymax = 0.1
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.05 0.01 0'
    radius = 0.01
    density = 2500
    gravity = '0 -9.81 0'
    normal_stiffness = 1e5
    normal_damping = 20
    wall_points = '0 0 0
                   0 0.0202 0'
    wall_normals = '0 1 0
                    0 -1 0'
    servo_walls = 1
    servo_forces = 1
    servo_gain = 0.02
    servo_max_velocity = 0.05
    skin = 0.001
    substeps = 500
    verify = true
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [top]
    type = KokkosWallForce
    cloud = cloud
    wall = 1
    component = y
  []
  [floor]
    type = KokkosWallForce
    cloud = cloud
    wall = 0
    component = y
  []
  [y]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = center_of_mass_y
  []
[]

[Executioner]
  type = Transient
  dt = 0.01
  num_steps = 150
[]

[Outputs]
  csv = true
[]
