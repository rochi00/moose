# The collision of restitution.i across a periodic face: the spheres move apart, each toward the
# nearer end of the periodic x range, and meet through their periodic images at x = +-0.2 with the
# same gap and closing speed, so the gold is that of restitution.i reflected about the face. The
# centers never cross the face (the overlap is a few percent of the radius), so this exercises
# contact through ghost images alone, across ranks when the spheres are on different ones.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 2
    xmin = -0.2
    xmax = 0.2
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '-0.14 0.5 0
                          0.14 0.5 0'
    initial_velocities = '-1 0 0
                           1 0 0'
    radius = 0.05
    density = 1000
    normal_stiffness = 1e5
    normal_damping = 64.7
    periodic = 'x'
    skin = 0.01
    substeps = 10000
    verify = true
    verify_neighbor_list = true
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
  dt = 0.005
  num_steps = 4
[]

[Outputs]
  csv = true
[]
