# Melting at a single temperature T_m with a solid of its own heat capacity and conductivity,
# against the two-phase Neumann similarity solution with different diffusivities (Carslaw and
# Jaeger, Conduction of Heat in Solids, 2nd ed., sec. 11.2). The solid has cp_s = 0.5 and
# k_s = 2, so its diffusivity is a_s = 4 against a_l = 1 of the liquid, and nu = sqrt(a_l / a_s)
# = 0.5. The front lies at x = 2 lambda sqrt(t), with lambda the root of
#   St_l / (exp(lambda^2) erf(lambda)) - St_s / (nu exp(lambda^2 nu^2) erfc(lambda nu))
#     = sqrt(pi) lambda
# for the Stefan numbers St_l = cp_l (T_w - T_m) / L = 1 and St_s = cp_s (T_m - T_i) / L = 0.25.
# The solid fluid properties use a reference temperature of 290 rather than the liquid's 300, so
# the enthalpies of the two objects differ by 5 at T_m; the mixture enthalpy aligns them there.

T_s = 300
T_l = ${T_s}
# Root of the equation above, to 15 digits
lambda = 0.49365853168693
nu = 0.5

!include stefan_base.i

[Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/heat]
  solid_fluid_properties = 'fp_s'
[]

[FluidProperties]
  [fp_s]
    type = TemperaturePressureFunctionFluidProperties
    rho = 1
    mu = 1
    k = 2
    cv = 0.5
    T_ref = 290
  []
[]

[Functions]
  [T_exact]
    type = ParsedFunction
    # t is floored so that the similarity variable stays finite at t = 0, where the profile is T_i
    expression = 'eta := x / (2 * sqrt(max(t, 1e-30)));
                  if(eta < lambda, T_w + (T_m - T_w) * erf(eta) / erf(lambda),
                                   T_i + (T_m - T_i) * (1 - erf(nu * eta)) / (1 - erf(nu * lambda)))'
    symbol_names = 'lambda nu T_w T_i T_m'
    symbol_values = '${lambda} ${nu} ${T_w} ${T_i} ${T_s}'
  []
  [melt_exact]
    type = ParsedFunction
    expression = '2 * lambda * sqrt(t)'
    symbol_names = 'lambda'
    symbol_values = '${lambda}'
  []
[]
