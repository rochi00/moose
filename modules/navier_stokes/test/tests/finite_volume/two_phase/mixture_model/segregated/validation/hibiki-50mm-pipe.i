##########################################################
# Vertical upward air-water bubbly flow in the 50.8 mm pipe of Hibiki, Ishii and Xiao (2001).
# The three axial stations of the experiment, z/D = 6.00, 30 and 54, are cut as internal
# sidesets; with 120 axial cells over 60 diameters they land exactly on cell faces.
# See hibiki-pipe-base.i for the model and the quantities reported.
##########################################################

# Geometry
D = 0.0508
R = 0.0254
L = 3.048       # 60 diameters
nx = 120

# Liquid superficial velocity of the first low void fraction case of the reference report
jf = 0.491

!include hibiki-pipe-base.i

[Mesh]
  # Cross-sections on which the area averages are taken
  [station_6]
    type = ParsedGenerateSideset
    input = gen
    combinatorial_geometry = 'abs(x - ${fparse 6 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD6'
  []
  [station_30]
    type = ParsedGenerateSideset
    input = station_6
    combinatorial_geometry = 'abs(x - ${fparse 30 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD30'
  []
  [station_54]
    type = ParsedGenerateSideset
    input = station_30
    combinatorial_geometry = 'abs(x - ${fparse 54 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD54'
  []
[]

[Postprocessors]
  # C0 = <alpha j> / (<alpha><j>), Vgj = <alpha (1-alpha) u_slip> / <alpha>
  [alpha_6]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD6'
  []
  [j_6]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD6'
  []
  [alpha_j_6]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD6'
  []
  [alpha_v_gj_6]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD6'
  []
  [C0_6]
    type = ParsedPostprocessor
    expression = 'alpha_j_6 / (alpha_6 * j_6)'
    pp_names = 'alpha_j_6 alpha_6 j_6'
  []
  [Vgj_6]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_6 / alpha_6'
    pp_names = 'alpha_v_gj_6 alpha_6'
  []
  [alpha_30]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD30'
  []
  [j_30]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD30'
  []
  [alpha_j_30]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD30'
  []
  [alpha_v_gj_30]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD30'
  []
  [C0_30]
    type = ParsedPostprocessor
    expression = 'alpha_j_30 / (alpha_30 * j_30)'
    pp_names = 'alpha_j_30 alpha_30 j_30'
  []
  [Vgj_30]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_30 / alpha_30'
    pp_names = 'alpha_v_gj_30 alpha_30'
  []
  [alpha_54]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD54'
  []
  [j_54]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD54'
  []
  [alpha_j_54]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD54'
  []
  [alpha_v_gj_54]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD54'
  []
  [C0_54]
    type = ParsedPostprocessor
    expression = 'alpha_j_54 / (alpha_54 * j_54)'
    pp_names = 'alpha_j_54 alpha_54 j_54'
  []
  [Vgj_54]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_54 / alpha_54'
    pp_names = 'alpha_v_gj_54 alpha_54'
  []
[]
