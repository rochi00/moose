# Melting at a single temperature T_m against the two-phase Neumann similarity solution
# (Carslaw and Jaeger, Conduction of Heat in Solids, 2nd ed., sec. 11.2). The front lies at
# x = 2 lambda sqrt(t), with lambda the root of
#   St_l / (exp(lambda^2) erf(lambda)) - St_s / (exp(lambda^2) erfc(lambda)) = sqrt(pi) lambda
# for the Stefan numbers St_l = cp (T_w - T_m) / L = 1 and St_s = cp (T_m - T_i) / L = 0.5.

T_s = 300
T_l = ${T_s}
# Root of the equation above, to 15 digits
lambda = 0.469850999707148

!include stefan_base.i

[Functions]
  [T_exact]
    type = ParsedFunction
    # t is floored so that the similarity variable stays finite at t = 0, where the profile is T_i
    expression = 'eta := x / (2 * sqrt(max(t, 1e-30)));
                  if(eta < lambda, T_w + (T_m - T_w) * erf(eta) / erf(lambda),
                                   T_i + (T_m - T_i) * (1 - erf(eta)) / (1 - erf(lambda)))'
    symbol_names = 'lambda T_w T_i T_m'
    symbol_values = '${lambda} ${T_w} ${T_i} ${T_s}'
  []
  [melt_exact]
    type = ParsedFunction
    expression = '2 * lambda * sqrt(t)'
    symbol_names = 'lambda'
    symbol_values = '${lambda}'
  []
[]
