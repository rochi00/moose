temporary_density = 1
density = 1
mass_flux_scale = 1

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 1
    dx = '0.5 0.5'
    ix = '1 1'
    subdomain_id = '1 2'
  []
[]

[Problem]
  linear_sys_names = 'temperature_system'
[]

[Variables]
  [T]
    type = MooseLinearVariableFVReal
    solver_sys = temperature_system
  []
[]

[ICs]
  [temperature]
    type = FunctionIC
    variable = T
    function = initial_temperature
  []
[]

[Functions]
  [initial_temperature]
    type = ParsedFunction
    expression = 'if(x < 0.5, 2, 4)'
  []
  [mass_flux]
    type = ParsedFunction
    expression = '${mass_flux_scale} * 4 * x * (1 - x)'
  []
[]

[FunctorMaterials]
  [enthalpy]
    type = ParsedFunctorMaterial
    property_name = specific_enthalpy
    functor_names = T
    expression = '4 * T + 1'
  []
  [constants]
    type = GenericFunctorMaterial
    prop_names = 'temporary_density density dh_dT'
    prop_values = '${temporary_density} ${density} 4'
  []
[]

[LinearFVKernels]
  [enthalpy_time]
    type = LinearFVLowMachEnthalpyTimeDerivative
    variable = T
    temporary_density = temporary_density
    density = density
    specific_enthalpy = specific_enthalpy
    dh_dT = dh_dT
  []
  [enthalpy_advection]
    type = LinearFVLowMachEnthalpyAdvection
    variable = T
    specific_enthalpy = specific_enthalpy
    dh_dT = dh_dT
    mass_flux_functor = mass_flux
  []
[]

[LinearFVBCs]
  [ends]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = T
    boundary = 'left right'
    functor = initial_temperature
  []
[]

[Postprocessors]
  [left_temperature]
    type = ElementAverageValue
    variable = T
    block = 1
  []
  [right_temperature]
    type = ElementAverageValue
    variable = T
    block = 2
  []
[]

[Executioner]
  type = Transient
  system_names = temperature_system
  scheme = bdf2
  num_steps = 1
  dt = 1
  l_tol = 1e-12
[]

[Outputs]
  csv = true
[]
