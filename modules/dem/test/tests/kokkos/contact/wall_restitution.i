# Oblique impact of one sphere on a fixed planar wall with the linear spring-dashpot model. The
# normal motion is the damped oscillator of restitution.i with the sphere's own mass as the
# effective mass, so e = 0.6385 and the contact lasts 7.3 ms (generate.py); the tangential velocity
# is untouched since the contact is frictionless. The wall at x = 0 lies inside the mesh so the
# sphere never leaves it.

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
    initial_positions = '-0.06 0.3 0'
    initial_velocity = '1 0.5 0'
    radius = 0.05
    density = 1000
    normal_stiffness = 1e5
    normal_damping = 64.7
    wall_points = '0 0 0'
    wall_normals = '-1 0 0'
    skin = 0.01
    substeps = 10000
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
  dt = 0.005
  num_steps = 4
[]

[Outputs]
  csv = true
[]
