[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
    xmin = 0
    xmax = 1
  []
[]

[Problem]
  linear_sys_names = 'mass_system'
[]

[Variables]
  [rho_adv]
    type = MooseLinearVariableFVReal
    solver_sys = mass_system
  []
[]

[ICs]
  [rho_adv]
    type = ConstantIC
    variable = rho_adv
    value = 1
  []
[]

[Functions]
  [volumetric_flux]
    type = ParsedFunction
    expression = 'if(x < 0.5, -1, 1)'
  []
[]

[FunctorMaterials]
  [eos]
    type = GenericFunctorMaterial
    prop_names = rho_eos
    prop_values = 1
  []
[]

[LinearFVKernels]
  [mass_time]
    type = LinearFVTimeDerivative
    variable = rho_adv
  []
  [mass_advection]
    type = LinearFVLowMachMassAdvection
    variable = rho_adv
    volumetric_face_flux = volumetric_flux
  []
[]

[LinearFVBCs]
  [inlet]
    type = LinearFVInletOutletScalarBC
    variable = rho_adv
    boundary = left
    face_flux = volumetric_flux
    backflow_value = 2
  []
  [outlet]
    type = LinearFVInletOutletScalarBC
    variable = rho_adv
    boundary = right
    face_flux = volumetric_flux
    backflow_value = 3
  []
[]

[UserObjects]
  [mass_flux]
    type = LowMachImplicitMassFlux
    system = mass_system
    mass_advection_kernel = mass_advection
    temporary_density = rho_adv
    eos_density = rho_eos
    mass_flux_name = low_mach_mass_flux
    maximum_density = 2
    execute_on = TIMESTEP_END
    execution_order_group = -1
  []
[]

[Postprocessors]
  [density]
    type = ElementAverageValue
    variable = rho_adv
  []
  [mass]
    type = ElementIntegralVariablePostprocessor
    variable = rho_adv
  []
  [inlet_mass_flux]
    type = SideIntegralFunctorPostprocessor
    boundary = left
    functor = low_mach_mass_flux
  []
  [outlet_mass_flux]
    type = SideIntegralFunctorPostprocessor
    boundary = right
    functor = low_mach_mass_flux
  []
[]

[Executioner]
  type = Transient
  system_names = mass_system
  scheme = bdf2
  num_steps = 1
  dt = 0.1
  l_tol = 1e-13
[]

[Outputs]
  csv = true
[]
