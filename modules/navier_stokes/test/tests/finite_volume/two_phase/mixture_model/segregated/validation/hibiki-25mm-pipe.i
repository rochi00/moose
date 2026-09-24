##########################################################
# Vertical upward air-water bubbly flow in the 25.4 mm pipe of Hibiki and Ishii (1999).
# The domain is 130 diameters long so that the z/D = 125 measuring station of Fig. 2 of Hibiki
# and Ishii (2003) sits inside the domain rather than on the outlet; the stations z/D = 6, 60 and
# 125 are cut as internal sidesets and land on cell faces with 260 axial cells.
# See hibiki-pipe-base.i for the model and the quantities reported.
##########################################################

# Geometry
D = 0.0254
R = 0.0127
L = 3.302       # 130 diameters
nx = 260

jf = 0.872

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
  [station_60]
    type = ParsedGenerateSideset
    input = station_6
    combinatorial_geometry = 'abs(x - ${fparse 60 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD60'
  []
  [station_125]
    type = ParsedGenerateSideset
    input = station_60
    combinatorial_geometry = 'abs(x - ${fparse 125 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD125'
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
  [alpha_60]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD60'
  []
  [j_60]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD60'
  []
  [alpha_j_60]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD60'
  []
  [alpha_v_gj_60]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD60'
  []
  [C0_60]
    type = ParsedPostprocessor
    expression = 'alpha_j_60 / (alpha_60 * j_60)'
    pp_names = 'alpha_j_60 alpha_60 j_60'
  []
  [Vgj_60]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_60 / alpha_60'
    pp_names = 'alpha_v_gj_60 alpha_60'
  []
  [alpha_125]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD125'
  []
  [j_125]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD125'
  []
  [alpha_j_125]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD125'
  []
  [alpha_v_gj_125]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD125'
  []
  [C0_125]
    type = ParsedPostprocessor
    expression = 'alpha_j_125 / (alpha_125 * j_125)'
    pp_names = 'alpha_j_125 alpha_125 j_125'
  []
  [Vgj_125]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_125 / alpha_125'
    pp_names = 'alpha_v_gj_125 alpha_125'
  []
[]
