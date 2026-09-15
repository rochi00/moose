# V3: head-on collision of two equal spheres with the linear spring-dashpot model. The relative
# motion is a damped linear oscillator, so the coefficient of restitution
#   e = exp(-zeta pi / sqrt(1 - zeta^2)),  zeta = gamma / (2 sqrt(k m_eff))
# and the contact duration pi / (omega_0 sqrt(1 - zeta^2)) are analytic; generate.py produces the
# gold state after the collision. With these parameters e = 0.5266 and the contact lasts 5.2 ms,
# resolved by 10400 substeps. The dashpot force is evaluated with the half-step velocity and the
# contact starts and ends on substep boundaries, both first order in dt, so the restitution error is
# 8e-5 here and 8e-6 at five times the substeps.

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
    initial_positions = '-0.06 0.5 0
                          0.06 0.5 0'
    initial_velocities = '1 0 0
                         -1 0 0'
    radius = 0.05
    density = 1000
    normal_stiffness = 1e5
    normal_damping = 64.7
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
