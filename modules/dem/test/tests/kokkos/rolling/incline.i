# V8: a sphere placed at rest on a plane inclined at 10 degrees (tan = 0.176) with friction and
# rolling resistance, started pressed in by the static normal overlap so the friction holds from
# the first substep. With friction = 0.5 and rolling_friction = 0.25 the friction holds the sphere
# and the rolling resistance torque balances the friction's, so it rests; with
# rolling_friction = 0.1 (incline_roll) the rolling spring saturates and the sphere rolls without
# slipping at g (sin - mu_r cos) / 1.4; with friction = 0.15 (incline_slide) the rolling resistance
# locks the rotation and the sphere slides at g (sin - mu cos) past the critical angle atan(mu).
# See generate.py. The dashpots damp the settling oscillations, of periods 2 ms and less, within
# 0.05 s.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 4
    xmax = 0.5
    ymin = -0.2
    ymax = 0.2
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.05173630609788627 0.1010307321708105 0'
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
    rolling_friction = 0.25
    wall_points = '0 0.1 0'
    wall_normals = '0.17364817766693033 0.984807753012208 0'
    skin = 0.001
    substeps = 2500
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
  num_steps = 10
[]

[Outputs]
  csv = true
  file_base = incline_rest_out
[]
