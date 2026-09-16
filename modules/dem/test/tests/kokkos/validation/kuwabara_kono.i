# Tier 1 validation: coefficient of restitution vs impact velocity with the Kuwabara-Kono
# viscoelastic damping of the Hertz model (dissipation_time = A). Two equal spheres collide
# head-on; the tests run closing speeds of 0.1, 0.3, 1, and 3 m/s (initial_velocities and the
# output name set from the test spec) and compare the rebound with the numerical integration of
# the viscoelastic collision (generate_kuwabara_kono.py): e falls from 0.809 to 0.661 across the
# decade and a half, following 1 - e ~ v^(1/5) at low speed, where the Tsuji damping would keep
# it constant.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 2
    xmin = -0.1
    xmax = 0.1
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '-0.01001 0.5 0
                          0.01001 0.5 0'
    initial_velocities = '0.5 0 0
                         -0.5 0 0'
    radius = 0.01
    density = 2500
    contact_model = hertz
    youngs_modulus = 1e8
    poissons_ratio = 0.3
    dissipation_time = 5e-5
    skin = 0.001
    substeps = 100
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
  dt = 2e-4
  num_steps = 12
[]

[Outputs]
  csv = true
[]
