# Tier 1 validation: Nesterenko's solitary wave in a chain of barely touching beads (Job, Melo,
# Sokolow, and Sen 2008, Eq. 2; Nesterenko 1983). Sixty 10 mm steel beads (Y = 203 GPa,
# nu = 0.3, 7780 kg/m3) lie along x, touching, with the elastic Hertz contact and nothing else;
# the first bead is set moving at 1 m/s, which launches a single solitary wave. The long
# wavelength theory gives the strain psi = psi_m cos^4((x - v t) / (R sqrt(10))) travelling at
# v = (6 / (5 pi rho theta))^(1/2) psi_m^(1/4), theta = 3 (1 - nu^2) / (4 Y), with the bead
# velocity v psi: the speed scales as the fifth root of the amplitude and the pulse spans about
# five bead diameters. nesterenko_check.py measures the speed, amplitude, and width from the
# state output every 2 microseconds and compares them with the theory.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 12
    ny = 1
    xmax = 0.6
    ymax = 0.1
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.0050000000000000001 0.05 0
                         0.014999999999999999 0.05 0
                         0.025000000000000001 0.05 0
                         0.035000000000000003 0.05 0
                         0.044999999999999998 0.05 0
                         0.055 0.05 0
                         0.065000000000000002 0.05 0
                         0.074999999999999997 0.05 0
                         0.085000000000000006 0.05 0
                         0.095000000000000001 0.05 0
                         0.105 0.05 0
                         0.115 0.05 0
                         0.125 0.05 0
                         0.13500000000000001 0.05 0
                         0.14499999999999999 0.05 0
                         0.155 0.05 0
                         0.16500000000000001 0.05 0
                         0.17500000000000002 0.05 0
                         0.185 0.05 0
                         0.19500000000000001 0.05 0
                         0.20500000000000002 0.05 0
                         0.215 0.05 0
                         0.22500000000000001 0.05 0
                         0.23500000000000001 0.05 0
                         0.245 0.05 0
                         0.255 0.05 0
                         0.26500000000000001 0.05 0
                         0.27500000000000002 0.05 0
                         0.28500000000000003 0.05 0
                         0.29499999999999998 0.05 0
                         0.30499999999999999 0.05 0
                         0.315 0.05 0
                         0.32500000000000001 0.05 0
                         0.33500000000000002 0.05 0
                         0.34500000000000003 0.05 0
                         0.35499999999999998 0.05 0
                         0.36499999999999999 0.05 0
                         0.375 0.05 0
                         0.38500000000000001 0.05 0
                         0.39500000000000002 0.05 0
                         0.40500000000000003 0.05 0
                         0.41500000000000004 0.05 0
                         0.42499999999999999 0.05 0
                         0.435 0.05 0
                         0.44500000000000001 0.05 0
                         0.45500000000000002 0.05 0
                         0.46500000000000002 0.05 0
                         0.47500000000000003 0.05 0
                         0.48499999999999999 0.05 0
                         0.495 0.05 0
                         0.505 0.05 0
                         0.51500000000000001 0.05 0
                         0.52500000000000002 0.05 0
                         0.53500000000000003 0.05 0
                         0.54500000000000004 0.05 0
                         0.55500000000000005 0.05 0
                         0.56500000000000006 0.05 0
                         0.57500000000000007 0.05 0
                         0.58499999999999996 0.05 0
                         0.59499999999999997 0.05 0'
    initial_velocities = '1 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0
                          0 0 0'
    radius = 0.005
    density = 7780
    contact_model = hertz
    youngs_modulus = 2.03e11
    poissons_ratio = 0.3
    restitution = 1
    skin = 0.0005
    substeps = 20
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
  dt = 2e-6
  num_steps = 400
[]

[Outputs]
  csv = true
[]
