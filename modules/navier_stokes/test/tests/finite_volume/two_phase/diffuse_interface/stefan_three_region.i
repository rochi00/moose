# Melting over a finite interval [T_s, T_l] with a linear liquid fraction, against the exact
# three-region (liquid, mushy, solid) similarity solution of E. N. Kondrashov, "Towards analytical
# solutions of the alloy solidification problem", arXiv:math-ph/0507046 (2006). With equal solid
# and liquid properties his eq. (11) makes the liquid fraction linear and the mushy diffusivity
# constant, km = k / (rho (cp + L / (T_l - T_s))) (eq. 17), so eqs. (27)-(29) give the profile in
# each region. The near-wall region is liquid here, as the wall is hot. The liquidus and solidus
# lie at x = 2 a sqrt(t) and x = 2 b sqrt(t), with a and b the roots of the interface conditions
# (34)-(35); as T_l - T_s -> 0 both tend to the Neumann constant of stefan_neumann.i at first order.

T_s = 300
T_l = 300.1
# Roots of Kondrashov's eqs. (34)-(35) for these data, to 15 digits
a = 0.419545440967726
b = 0.493031972947813

!include stefan_base.i

km = ${fparse 1 / (1 + L / (T_l - T_s))}

[Functions]
  [T_exact]
    type = ParsedFunction
    # t is floored so that the similarity variable stays finite at t = 0, where the profile is T_i
    expression = 'eta := x / (2 * sqrt(max(t, 1e-30)));
                  if(eta < a, T_w + (T_l - T_w) * erf(eta) / erf(a),
                  if(eta < b, T_s + (T_l - T_s) * (erf(b / sqrt(km)) - erf(eta / sqrt(km)))
                                                / (erf(b / sqrt(km)) - erf(a / sqrt(km))),
                              T_i + (T_s - T_i) * (1 - erf(eta)) / (1 - erf(b))))'
    symbol_names = 'a b km T_w T_i T_s T_l'
    symbol_values = '${a} ${b} ${km} ${T_w} ${T_i} ${T_s} ${T_l}'
  []
  # The liquid depth plus the integral of the linear liquid fraction over the mushy region,
  # using the antiderivative s erf(s) + exp(-s^2) / sqrt(pi) of erf
  [melt_exact]
    type = ParsedFunction
    expression = 's1 := a / sqrt(km);
                  s2 := b / sqrt(km);
                  F1 := s1 * erf(s1) + exp(-s1 * s1) / sqrt(pi);
                  F2 := s2 * erf(s2) + exp(-s2 * s2) / sqrt(pi);
                  2 * a * sqrt(t)
                  + 2 * sqrt(km * t) * (erf(s2) * (s2 - s1) - (F2 - F1)) / (erf(s2) - erf(s1))'
    symbol_names = 'a b km'
    symbol_values = '${a} ${b} ${km}'
  []
[]
