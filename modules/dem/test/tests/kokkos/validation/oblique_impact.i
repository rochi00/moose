# Tier 1 validation: oblique impact of a sphere on a plane (the setting of Kharaz, Gorham, and
# Salman 2001: a 5 mm aluminium oxide sphere on a soda-lime glass anvil at 3.85 m/s, e_n = 0.98,
# mu = 0.092). The single material of the cloud is chosen to give the effective moduli of the
# alumina-glass pair: E* = 63 GPa and G* = 13.5 GPa. The impact angle is set from the initial
# velocity by the sweep (oblique_impact_sweep.py) or the test spec; the rebound spin and the
# tangential velocity are read from the state after the 13 microsecond contact.
#
# Above the critical angle, tan(theta) > (7/2) mu (1 + e_n) = 0.64 (theta > 32.5 degrees), the
# contact slides throughout and the rebound is analytic: v_t' = v_t - mu (1 + e_n) v_n and
# omega' = (5/2) mu (1 + e_n) v_n / r, whatever the tangential model. Below it the contact sticks
# for part of the impact and the rebound depends on the tangential model (Maw, Barber, and
# Fawcett 1976).

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
    xmin = -0.01
    xmax = 0.01
    ymax = 0.02
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0 0.002501 0'
    initial_velocity = '3.334 -1.925 0'
    radius = 0.0025
    density = 3900
    contact_model = hertz
    youngs_modulus = 1.187e11
    poissons_ratio = 0.24
    restitution = 0.98
    friction = 0.092
    wall_points = '0 0 0'
    wall_normals = '0 1 0'
    skin = 0.0002
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
  dt = 1e-5
  num_steps = 3
[]

[Outputs]
  csv = true
[]
