##########################################################
# Common part of the three-dimensional T-junction cases: a square main channel with a side branch,
# Re = 100 on the inlet velocity and channel width, carrying a dispersed phase lighter than the
# continuous one so that gravity separates the two and the drift flux matters. The including input
# picks the discretization: the linear FV segregated solver (simple.i, pimple.i) or the nonlinear
# FV Newton solver (newton.i).
##########################################################

width = 0.5
half_width = '${fparse 0.5 * width}'
main_length = 2
branch_length = 1
junction_end = '${fparse main_length + width}'
branch_end = '${fparse branch_length + half_width}'

# A 0.01 m/s inlet and Re = 100 give mu = rho * U * D_h / Re = 0.05 Pa s. These are the
# continuous phase properties; the mixture density and viscosity the momentum equation carries
# are built by the two-phase Physics from them and from the dispersed phase below.
rho = 1000
inlet_velocity = 0.01
mu = 0.05

# Dispersed phase: ten times lighter, carried in at a 10% volume fraction as 1 mm particles.
# Gravity acts across the channel (-z), so the light phase rises towards the top wall as it is
# carried along. The Stokes rise velocity of a 1 mm particle at this density difference is
# 0.0098 m/s, the same order as the flow, so the phase stratifies within one channel length and
# the drift flux is a leading term rather than a correction.
#
# Both solvers solve for the dynamic pressure, with the hydrostatic head of the inlet mixture
# absorbed into the pressure, so the outlets can stay at a uniform zero pressure. Gravity then
# acts on the slip closure (the separation) but the buoyancy of the mixture density variation is
# not carried in the momentum equation; that is the same simplification on both sides.
rho_d = 100
mu_d = 0.005
dp = 0.001
inlet_phase_2 = 0.1
gravity = '0 0 -9.81'

# Required by the mixture Physics, unused without an energy equation
cp = 1
k = 1
cp_d = 1
k_d = 1

[Mesh]
  [main]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = ${main_length}
    ymin = -${half_width}
    ymax = ${half_width}
    zmin = -${half_width}
    zmax = ${half_width}
    # Four cells across the channel give a uniform 0.125 m medium-coarse mesh.
    nx = 16
    ny = 4
    nz = 4
  []
  [name_main_boundaries]
    type = RenameBoundaryGenerator
    input = main
    old_boundary = 'left right bottom top back front'
    new_boundary = 'main_left main_right main_bottom main_top main_back main_front'
  []

  [junction]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = ${main_length}
    xmax = ${junction_end}
    ymin = -${half_width}
    ymax = ${half_width}
    zmin = -${half_width}
    zmax = ${half_width}
    nx = 4
    ny = 4
    nz = 4
  []
  [name_junction_boundaries]
    type = RenameBoundaryGenerator
    input = junction
    old_boundary = 'left right bottom top back front'
    new_boundary = 'junction_left junction_right junction_bottom junction_top junction_back junction_front'
  []

  [minus_y_branch]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = ${main_length}
    xmax = ${junction_end}
    ymin = -${branch_end}
    ymax = -${half_width}
    zmin = -${half_width}
    zmax = ${half_width}
    nx = 4
    ny = 8
    nz = 4
  []
  [name_minus_y_boundaries]
    type = RenameBoundaryGenerator
    input = minus_y_branch
    old_boundary = 'left right bottom top back front'
    new_boundary = 'minus_y_left minus_y_right minus_y_bottom minus_y_top minus_y_back minus_y_front'
  []

  [plus_y_branch]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = ${main_length}
    xmax = ${junction_end}
    ymin = ${half_width}
    ymax = ${branch_end}
    zmin = -${half_width}
    zmax = ${half_width}
    nx = 4
    ny = 8
    nz = 4
  []
  [name_plus_y_boundaries]
    type = RenameBoundaryGenerator
    input = plus_y_branch
    old_boundary = 'left right bottom top back front'
    new_boundary = 'plus_y_left plus_y_right plus_y_bottom plus_y_top plus_y_back plus_y_front'
  []

  [stitch]
    type = StitchMeshGenerator
    inputs = 'name_main_boundaries name_junction_boundaries name_minus_y_boundaries name_plus_y_boundaries'
    stitch_boundaries_pairs = 'main_right junction_left;
                               junction_bottom minus_y_top;
                               junction_top plus_y_bottom'
    clear_stitched_boundary_ids = true
  []

  [name_external_boundaries]
    type = RenameBoundaryGenerator
    input = stitch
    old_boundary = 'main_left minus_y_bottom plus_y_top
                    main_bottom main_top main_back main_front
                    junction_right junction_back junction_front
                    minus_y_left minus_y_right minus_y_back minus_y_front
                    plus_y_left plus_y_right plus_y_back plus_y_front'
    new_boundary = 'inlet outlet_minus_y outlet_plus_y
                    walls walls walls walls
                    walls walls walls
                    walls walls walls walls
                    walls walls walls walls'
  []
[]

[Postprocessors]
  # Flow split and phase split between the two branches, and the pressure drop
  [vy_minus]
    type = SideAverageFunctorPostprocessor
    functor = vel_y
    boundary = outlet_minus_y
  []
  [vy_plus]
    type = SideAverageFunctorPostprocessor
    functor = vel_y
    boundary = outlet_plus_y
  []
  [alpha_minus]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = outlet_minus_y
  []
  [alpha_plus]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = outlet_plus_y
  []
  [p_in]
    type = SideAverageFunctorPostprocessor
    functor = pressure
    boundary = inlet
  []
  [max_slip_x]
    type = ElementExtremeFunctorValue
    functor = vel_slip_x
  []
  [max_alpha]
    type = ElementExtremeFunctorValue
    functor = phase_2
  []
  [min_alpha]
    type = ElementExtremeFunctorValue
    functor = phase_2
    value_type = min
  []
[]

[Outputs]
  # The fields at the converged steady state only; the csv carries the transient
  [exodus]
    type = Exodus
    execute_on = FINAL
  []
  csv = true
[]
