# Oblique collision of two equal spheres with the Hertz-Mindlin model and Coulomb friction: the
# spheres approach head-on with an offset of 0.4 diameters, so the contact has a tangential slip
# that the friction torque converts into equal spins on both (the tangential force is equal and
# opposite and the torque arms are equal). The gold is the serial result; the test checks that the
# contact history advanced on both ranks of a pair across a partition, from the forwarded ghost
# positions, velocities, and spins, gives the same collision on any number of ranks.

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
    initial_positions = '-0.06 0.48 0
                          0.06 0.52 0'
    initial_velocities = '1 0 0
                         -1 0 0'
    radius = 0.05
    density = 1000
    contact_model = hertz
    youngs_modulus = 1e7
    poissons_ratio = 0.3
    restitution = 0.8
    friction = 0.3
    skin = 0.01
    substeps = 500
    verify = true
    verify_neighbor_list = true
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
  dt = 0.005
  num_steps = 6
[]

[Outputs]
  csv = true
[]
