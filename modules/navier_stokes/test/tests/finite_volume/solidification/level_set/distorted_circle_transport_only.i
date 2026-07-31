# Forward-reverse vortex transport without redistancing or AMR. This isolates the temporal and
# spatial accuracy of LinearFVMaterialAdvection from interface correction and mesh transfer.

cells = 32
finest_h = ${fparse 1.0 / cells}
interface_width = ${fparse 1.5 / cells}
radius = 0.15
center_x = 0.5
center_y = 0.75
half_time = 0.625
end_time = 1.25
velocity_scale = 1

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${cells}
    ny = ${cells}
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
[]

[Problem]
  linear_sys_names = transport_system
[]

[Variables]
  [level_set_phi]
    type = MooseLinearVariableFVReal
    solver_sys = transport_system
  []
[]

[Functions]
  [initial_circle]
    type = ParsedFunction
    expression = '${radius} - sqrt((x - ${center_x})^2 + (y - ${center_y})^2)'
  []
  [vortex_u]
    type = ParsedFunction
    expression = 's:=if(t <= ${half_time}, 1, -1); '
                 '-${velocity_scale} * s * sin(pi*x)^2 * sin(2*pi*y)'
  []
  [vortex_v]
    type = ParsedFunction
    expression = 's:=if(t <= ${half_time}, 1, -1); '
                 '${velocity_scale} * s * sin(pi*y)^2 * sin(2*pi*x)'
  []
[]

[ICs]
  [circle]
    type = FunctionIC
    variable = level_set_phi
    function = initial_circle
  []
[]

[FunctorMaterials]
  [masked_level_set]
    type = ParsedFunctorMaterial
    property_name = masked_level_set
    functor_names = level_set_phi
    functor_symbols = phi
    expression = 'd:=${radius} - sqrt((x - ${center_x})^2 + (y - ${center_y})^2); '
                 'if(abs(d) <= 3 * ${finest_h}, phi, d)'
  []
  [band_level_set_absolute_error]
    type = ParsedFunctorMaterial
    property_name = band_level_set_absolute_error
    functor_names = masked_level_set
    functor_symbols = phi
    expression = 'd:=${radius} - sqrt((x - ${center_x})^2 + (y - ${center_y})^2); '
                 'abs(phi - d)'
  []
[]

[LinearFVKernels]
  [time]
    type = LinearFVTimeDerivative
    variable = level_set_phi
  []
  [advection]
    type = LinearFVMaterialAdvection
    variable = level_set_phi
    u = vortex_u
    v = vortex_v
    advected_interp_method = venkatakrishnan
  []
[]

[Postprocessors]
  [band_level_set_l2_error]
    type = ElementL2FunctorError
    approximate = masked_level_set
    exact = initial_circle
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [band_level_set_linf_error]
    type = ElementExtremeFunctorValue
    functor = band_level_set_absolute_error
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  system_names = transport_system
  dt = ${finest_h}
  end_time = ${end_time}
  l_tol = 1e-12
  [TimeIntegrator]
    type = BDF2
    restart_times = '${half_time}'
  []
[]

[Outputs]
  csv = true
  exodus = false
  file_base = distorted_circle_transport_only
  execute_on = 'INITIAL TIMESTEP_END'
  sync_times = '${half_time}'
  [console]
    type = Console
    execute_postprocessors_on = 'INITIAL FINAL'
  []
[]
