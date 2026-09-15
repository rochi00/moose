# V7: a sphere resting on a wall under gravity with the Hertz-Mindlin model and partial slip is
# given a small tangential velocity. The normal force stays m g, so the tangential force follows
# the Mindlin-Deresiewicz loading curve mu P [1 - (1 - delta / delta_max)^(3/2)] of the contact
# point displacement delta, and the sphere oscillates on it (elastic on unloading) with a period
# of 5 ms, reaching 0.7 delta_max: well into partial slip, short of gross sliding. The gold is
# the reference integration of generate_mindlin.py; the forces are output to compare the curve
# itself. restitution = 1 gives no damping, so the sphere starts at the static overlap.

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
    initial_positions = '0.05 0.009973008924843728 0'
    initial_velocity = '0.008047460718913348 0 0'
    radius = 0.01
    density = 2500
    gravity = '0 -9.81 0'
    contact_model = hertz
    youngs_modulus = 1e7
    poissons_ratio = 0.3
    restitution = 1
    friction = 0.3
    partial_slip = true
    wall_points = '0 0 0'
    wall_normals = '0 1 0'
    skin = 0.001
    substeps = 10000
    verify = true
    verify_pair_states = true
    execute_on = TIMESTEP_END
  []
[]

[VectorPostprocessors]
  [state]
    type = KokkosParticleState
    cloud = cloud
    output_forces = true
    execute_on = TIMESTEP_END
  []
[]

[Executioner]
  type = Transient
  dt = 0.002
  num_steps = 10
[]

[Outputs]
  csv = true
[]
