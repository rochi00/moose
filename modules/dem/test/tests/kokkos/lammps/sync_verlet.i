# LAMMPS examples/granular/in.sync_verlet: a 0.5 mm sphere (2500 kg/m3) dropped onto two frozen
# 10 mm spheres, pair granular hooke 1e4 with restitution 0.5 (damping coeff_restitution: the
# dashpot 2 beta sqrt(k m_eff) = 0.017431 N s/m with beta = -ln e / sqrt(pi^2 + ln^2 e) and the
# fine sphere's mass, its frozen partners being immovable), linear_history 8235 with no
# tangential damping, mu = 0.5, limit_damping; walls at z = +-0.02. With pair granular's
# conventions (springs rescaled on rotation, torque arm at the contact plane) the plain-Verlet
# LAMMPS run (lammps/in.sync_verlet with variant plain) is reproduced to roundoff. The
# original synchronized_verlet variant differs in its second history projection and, this
# trajectory being sensitive (the sphere may or may not escape the crevice), parts from the
# plain one after 0.15 s; LAMMPS's own shipped log parts from a current LAMMPS run there too
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 12
    ny = 8
    nz = 8
    xmin = -0.03
    xmax = 0.03
    ymin = -0.02
    ymax = 0.02
    zmin = -0.02
    zmax = 0.02
  []
[]
[Problem]
  solve = false
[]
[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.005 2.7937481665024794e-05 0.0016005372526493155
                         0.01 0 0
                         0 0 0'
    initial_radii = '0.00025 0.005 0.005'
    initial_frozen = '1 2'
    radius = 0.005
    density = 2500
    gravity = '0 0 -9.81'
    normal_stiffness = 1e4
    normal_damping = 0.017430441966303
    tangential_stiffness = 8235
    tangential_damping = 0
    friction = 0.5
    limit_damping = true
    deformed_torque_arm = true
    wall_points = '0 0 -0.02
                   0 0 0.02'
    wall_normals = '0 0 1
                    0 0 -1'
    periodic = 'x y'
    skin = 0.0005
    substeps = 60000
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
  dt = 0.06
  num_steps = 5
[]
[Outputs]
  csv = true
[]
