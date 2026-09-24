##########################################################
# Radial void fraction and axial velocity profiles in the 50.8 mm pipe of Hibiki, Ishii and Xiao
# (2001), at the z/D = 6, 30 and 53.5 stations. This is the quantity Swiderski et al. validate against
# for this facility, and the one a flat void profile cannot reproduce.
# See hibiki-pipe-base.i for the model.
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
  [station_53p5]
    type = ParsedGenerateSideset
    input = station_30
    combinatorial_geometry = 'abs(x - ${fparse 53.5 * D}) < 1e-8'
    normal = '1 0 0'
    new_sideset_name = 'zD53p5'
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
  [alpha_53p5]
    type = SideAverageFunctorPostprocessor
    functor = phase_2
    boundary = 'zD53p5'
  []
  [j_53p5]
    type = SideAverageFunctorPostprocessor
    functor = j_axial
    boundary = 'zD53p5'
  []
  [alpha_j_53p5]
    type = SideAverageFunctorPostprocessor
    functor = alpha_j
    boundary = 'zD53p5'
  []
  [alpha_v_gj_53p5]
    type = SideAverageFunctorPostprocessor
    functor = alpha_v_gj
    boundary = 'zD53p5'
  []
  [C0_53p5]
    type = ParsedPostprocessor
    expression = 'alpha_j_53p5 / (alpha_53p5 * j_53p5)'
    pp_names = 'alpha_j_53p5 alpha_53p5 j_53p5'
  []
  [Vgj_53p5]
    type = ParsedPostprocessor
    expression = 'alpha_v_gj_53p5 / alpha_53p5'
    pp_names = 'alpha_v_gj_53p5 alpha_53p5'
  []
[]

# Radial profiles at the two experimental stations. This is the quantity Swiderski et al. validate
# against for this same facility, and the one a flat void profile cannot reproduce.
[Positions]
  # The 16 radial cell centroids at each experimental station
  [radial_zD6]
    type = InputPositions
    positions = '0.304800 0.000794 0' '0.304800 0.002381 0' '0.304800 0.003969 0' '0.304800 0.005556 0' '0.304800 0.007144 0' '0.304800 0.008731 0' '0.304800 0.010319 0' '0.304800 0.011906 0' '0.304800 0.013494 0' '0.304800 0.015081 0' '0.304800 0.016669 0' '0.304800 0.018256 0' '0.304800 0.019844 0' '0.304800 0.021431 0' '0.304800 0.023019 0' '0.304800 0.024606 0'
  []
  [radial_zD53p5]
    type = InputPositions
    positions = '2.717800 0.000794 0' '2.717800 0.002381 0' '2.717800 0.003969 0' '2.717800 0.005556 0' '2.717800 0.007144 0' '2.717800 0.008731 0' '2.717800 0.010319 0' '2.717800 0.011906 0' '2.717800 0.013494 0' '2.717800 0.015081 0' '2.717800 0.016669 0' '2.717800 0.018256 0' '2.717800 0.019844 0' '2.717800 0.021431 0' '2.717800 0.023019 0' '2.717800 0.024606 0'
  []
[]

[VectorPostprocessors]
  # Radial void fraction and axial velocity, the quantity Swiderski et al. validate against for
  # this facility. Sampled as functors because the variables live in a linear FV system.
  [prof_zD6]
    type = PositionsFunctorValueSampler
    positions = radial_zD6
    functors = 'phase_2 vel_x'
    discontinuous = true
    sort_by = y
    execute_on = 'FINAL'
  []
  [prof_zD53p5]
    type = PositionsFunctorValueSampler
    positions = radial_zD53p5
    functors = 'phase_2 vel_x'
    discontinuous = true
    sort_by = y
    execute_on = 'FINAL'
  []
[]
