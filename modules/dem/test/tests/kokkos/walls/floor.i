# V13: a sphere dropped obliquely on the floor of a 2D mesh, the floor being the 'bottom' sideset
# turned into a wall (segments in 2D), with the linear spring-dashpot model and friction. The
# floor_analytic test runs the same drop against the analytic plane y = 0 instead; the two must
# agree to roundoff, the sideset wall reducing to the plane. The sphere bounces twice in the run.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 4
    xmax = 0.4
    ymax = 0.2
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.05 0.08 0'
    initial_velocity = '0.5 -0.3 0'
    radius = 0.01
    density = 2500
    gravity = '0 -9.81 0'
    normal_stiffness = 1e5
    normal_damping = 1
    tangential_stiffness = 1e5
    tangential_damping = 1
    friction = 0.3
    wall_boundaries = 'bottom'
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
  dt = 0.05
  num_steps = 6
[]

[Outputs]
  csv = true
[]
