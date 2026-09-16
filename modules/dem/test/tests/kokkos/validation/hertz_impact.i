# Tier 1 validation: normal elastic impact of two equal spheres with the Hertz model at a closing
# speed of 1 m/s. The force history, the maximum overlap, and the contact duration follow Hertz's
# impact solution (generate_hertz_impact.py): d_max = 2.76e-4 m, t_c = 0.832 ms, F_max = 23.7 N.
# The state is output every 1e-5 s, some 80 points over the contact, and compared at eight of them
# against the reference integration; five substeps per output resolve the contact with 400.

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
    restitution = 1
    skin = 0.001
    substeps = 5
    verify = true
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
  dt = 1e-5
  num_steps = 90
[]

[Outputs]
  csv = true
[]
