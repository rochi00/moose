# V12, wedge: an L-shaped mesh (a 2x2 arrangement of unit blocks with the top-right one deleted)
# whose new boundary, the two segments meeting at the reentrant corner (1, 1), is a wall. A
# sphere in the bottom-left block hits the corner head-on along the diagonal: the two segments
# report the same closest point, the vertex, reduced to one contact along the line to the
# center, so the collision is the sphere-wall damped oscillator with e = 0.6385 and the sphere
# comes straight back (generate.py). A doubled contact would give a stiffer, different bounce.

[Mesh]
  [blocks]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1 1'
    dy = '1 1'
    ix = '4 4'
    iy = '4 4'
    subdomain_id = '0 1 2 3'
  []
  [delete]
    type = BlockDeletionGenerator
    input = blocks
    block = 3
    new_boundary = corner
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '0.8 0.8 0'
    initial_velocity = '0.70710678118654746 0.70710678118654746 0'
    radius = 0.05
    density = 1000
    normal_stiffness = 1e5
    normal_damping = 64.7
    wall_boundaries = 'corner'
    skin = 0.01
    substeps = 20000
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
  dt = 0.05
  num_steps = 8
[]

[Outputs]
  csv = true
[]
