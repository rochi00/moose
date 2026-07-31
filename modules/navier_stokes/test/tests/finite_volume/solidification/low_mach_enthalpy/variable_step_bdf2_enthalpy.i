!include conservative_enthalpy.i

[LinearFVKernels]
  active = 'enthalpy_time reaction'
  [reaction]
    type = LinearFVReaction
    variable = T
  []
[]

[Executioner]
  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '0 0.1 0.3 0.35'
  []
[]
