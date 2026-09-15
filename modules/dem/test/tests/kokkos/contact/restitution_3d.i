# restitution.i along the (1, 1, 1) diagonal of a three-dimensional HEX8 box centered on the
# origin: the analytic solution is the same and every position and velocity component is
# exercised, as is contact across 3D partitions. See generate.py for the gold.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 4
    nz = 4
    xmin = -0.2
    xmax = 0.2
    ymin = -0.2
    ymax = 0.2
    zmin = -0.2
    zmax = 0.2
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '-0.034641016151377546 -0.034641016151377546 -0.034641016151377546
                          0.034641016151377546 0.034641016151377546 0.034641016151377546'
    initial_velocities = '0.57735026918962584 0.57735026918962584 0.57735026918962584
                          -0.57735026918962584 -0.57735026918962584 -0.57735026918962584'
    radius = 0.05
    density = 1000
    normal_stiffness = 1e5
    normal_damping = 64.7
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
