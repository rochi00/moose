# V5: elastic granular gas of 64 spheres in a fully periodic unit box, from a lattice with
# pseudo-random velocities of zero total momentum (generate.py). Without damping the total energy
# (kinetic plus spring) is conserved to the integration error and the linear and angular momenta
# exactly, on any number of ranks, while particles wrap through the periodic faces and collide
# across them through ghost images. The contact lasts 2.4 ms, resolved by 95 substeps. The
# coordination number, contact count, virial pressure, and load imbalance are reported as
# diagnostics but not compared: they fluctuate with the collisions and, at this size, with the
# roundoff of the summation order.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 5
    ny = 5
    nz = 5
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.125 0.125 0.125
                         0.125 0.125 0.375
                         0.125 0.125 0.625
                         0.125 0.125 0.875
                         0.125 0.375 0.125
                         0.125 0.375 0.375
                         0.125 0.375 0.625
                         0.125 0.375 0.875
                         0.125 0.625 0.125
                         0.125 0.625 0.375
                         0.125 0.625 0.625
                         0.125 0.625 0.875
                         0.125 0.875 0.125
                         0.125 0.875 0.375
                         0.125 0.875 0.625
                         0.125 0.875 0.875
                         0.375 0.125 0.125
                         0.375 0.125 0.375
                         0.375 0.125 0.625
                         0.375 0.125 0.875
                         0.375 0.375 0.125
                         0.375 0.375 0.375
                         0.375 0.375 0.625
                         0.375 0.375 0.875
                         0.375 0.625 0.125
                         0.375 0.625 0.375
                         0.375 0.625 0.625
                         0.375 0.625 0.875
                         0.375 0.875 0.125
                         0.375 0.875 0.375
                         0.375 0.875 0.625
                         0.375 0.875 0.875
                         0.625 0.125 0.125
                         0.625 0.125 0.375
                         0.625 0.125 0.625
                         0.625 0.125 0.875
                         0.625 0.375 0.125
                         0.625 0.375 0.375
                         0.625 0.375 0.625
                         0.625 0.375 0.875
                         0.625 0.625 0.125
                         0.625 0.625 0.375
                         0.625 0.625 0.625
                         0.625 0.625 0.875
                         0.625 0.875 0.125
                         0.625 0.875 0.375
                         0.625 0.875 0.625
                         0.625 0.875 0.875
                         0.875 0.125 0.125
                         0.875 0.125 0.375
                         0.875 0.125 0.625
                         0.875 0.125 0.875
                         0.875 0.375 0.125
                         0.875 0.375 0.375
                         0.875 0.375 0.625
                         0.875 0.375 0.875
                         0.875 0.625 0.125
                         0.875 0.625 0.375
                         0.875 0.625 0.625
                         0.875 0.625 0.875
                         0.875 0.875 0.125
                         0.875 0.875 0.375
                         0.875 0.875 0.625
                         0.875 0.875 0.875'
    initial_velocities = '0.47981 -0.632034 0.742705
                          -1.713732 0.215007 0.001528
                          0.269083 -0.371472 -0.930469
                          -0.644077 1.451049 -1.266258
                          -0.949558 0.722832 1.201482
                          1.810439 1.350992 -0.100109
                          0.015189 0.650681 -0.95717
                          0.674746 1.013843 1.960933
                          -0.824719 -0.069461 0.876469
                          0.818773 -1.160178 -1.894316
                          0.989504 -1.686466 0.432593
                          -1.158017 0.374555 0.103739
                          -0.551004 -1.120825 0.626691
                          0.725371 -0.631522 1.928777
                          1.197324 -0.290923 0.862356
                          -1.638741 0.57977 0.239167
                          0.658728 1.764657 -0.369641
                          1.188447 0.62257 -1.892524
                          -0.640874 -1.413017 0.283873
                          -0.660071 -1.251812 1.257668
                          -1.852136 -1.508389 -1.169335
                          1.452608 0.345418 -0.001535
                          -0.635821 -1.451541 -0.143492
                          -0.906918 -0.34295 -1.466451
                          1.640739 1.316441 -0.219975
                          0.12083 -1.766932 1.883934
                          -1.502406 -0.455223 1.037571
                          0.739731 1.476995 -1.224929
                          0.653286 -1.642421 2.007753
                          -2.018883 -1.209748 2.020896
                          0.565817 -1.694499 -0.99142
                          0.416191 -1.219633 1.679236
                          -0.385386 1.59793 0.036176
                          -1.487151 1.272282 -1.942309
                          1.045549 -0.678608 -1.437524
                          -1.268351 -0.052048 -0.983744
                          0.824531 1.004499 -0.82982
                          -1.302184 0.837881 0.442679
                          1.381112 -0.992989 -1.423237
                          0.090914 2.099151 -0.231286
                          -0.628004 -1.609793 -0.814796
                          0.760943 -0.962044 -0.549885
                          -0.444749 -0.331998 1.197592
                          1.509363 0.670597 -0.049586
                          0.139518 0.642897 1.038663
                          1.465697 -0.2347 -0.551782
                          -1.94163 -0.101938 -0.14022
                          -1.454936 0.768202 1.729915
                          1.018489 1.946522 1.72981
                          1.589188 0.790239 0.034425
                          0.381209 2.030606 -1.732235
                          -2.026903 -0.582351 0.483751
                          0.460447 -1.380473 0.433041
                          0.20195 1.007215 1.662832
                          1.155956 -0.729889 -1.537495
                          0.827057 1.838316 0.404345
                          1.48929 -1.804583 1.029629
                          0.377419 1.517304 -0.611178
                          1.0708 -0.871847 -1.900696
                          1.432742 0.865639 1.287606
                          -1.700583 -1.629497 -1.837099
                          -1.043911 0.671674 -1.644116
                          -0.317083 1.789129 -0.938412
                          -1.9409619999999999 -1.3530890000000007 1.125209'
    radius = 0.03
    density = 1000
    normal_stiffness = 1e5
    periodic = 'x y z'
    substeps = 400
    verify = true
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [num_particles]
    type = KokkosParticleCloudValue
    cloud = cloud
    execute_on = 'INITIAL TIMESTEP_END'
    value = num_particles
  []
  [total_energy]
    type = KokkosParticleCloudValue
    cloud = cloud
    execute_on = 'INITIAL TIMESTEP_END'
    value = total_energy
  []
  [linear_momentum]
    type = KokkosParticleCloudValue
    cloud = cloud
    execute_on = 'INITIAL TIMESTEP_END'
    value = linear_momentum
  []
  [angular_momentum]
    type = KokkosParticleCloudValue
    cloud = cloud
    execute_on = 'INITIAL TIMESTEP_END'
    value = angular_momentum
  []
  [num_contacts]
    type = KokkosParticleCloudValue
    cloud = cloud
    execute_on = 'INITIAL TIMESTEP_END'
    value = num_contacts
  []
  [coordination_number]
    type = KokkosParticleCloudValue
    cloud = cloud
    execute_on = 'INITIAL TIMESTEP_END'
    value = coordination_number
  []
  [virial_pressure]
    type = KokkosParticleCloudValue
    cloud = cloud
    execute_on = 'INITIAL TIMESTEP_END'
    value = virial_pressure
  []
  [load_imbalance]
    type = KokkosParticleCloudValue
    cloud = cloud
    execute_on = 'INITIAL TIMESTEP_END'
    value = load_imbalance
  []
[]

[Executioner]
  type = Transient
  dt = 0.01
  num_steps = 30
[]

[Outputs]
  csv = true
[]
