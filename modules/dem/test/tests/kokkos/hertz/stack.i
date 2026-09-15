# V4: two equal spheres stacked on a wall under gravity with the Hertz contact model. Starting
# exactly touching, they settle to the overlaps at which (4/3) E* sqrt(r_eff) delta^(3/2) balances
# the weight above each contact: 2 m g against the wall (r_eff = r) and m g between the spheres
# (r_eff = r / 2), 43 and 34 micrometers here (generate.py). The contacts oscillate with periods
# of 7.6 and 6.8 ms, resolved by 400 substeps each, and the damping set by restitution = 0.05
# settles the positions to the gold within 0.2 s.

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
    initial_positions = '0.05 0.01 0
                         0.05 0.03 0'
    radius = 0.01
    density = 2500
    gravity = '0 -9.81 0'
    contact_model = hertz
    youngs_modulus = 1e7
    poissons_ratio = 0.3
    restitution = 0.05
    wall_points = '0 0 0'
    wall_normals = '0 1 0'
    skin = 0.001
    substeps = 1000
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
  dt = 0.02
  num_steps = 10
[]

[Outputs]
  csv = true
[]
