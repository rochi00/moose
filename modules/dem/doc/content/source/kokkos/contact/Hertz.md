# Hertz

Hertz normal contact model with the damping of Tsuji, Tanaka, and Ishida (Powder Technology 71,
1992) used by [KokkosParticleCloud.md]
when `contact_model = hertz` (plan layer L4). Like [LinearSpringDashpot.md] it is a compile-time
policy struct whose force the contact kernel inlines (plan decision D5). The normal force on
particle $i$ from particle $j$ is

!equation
F_n = \frac{4}{3} E^* \sqrt{r_\text{eff}} \, \delta^{3/2}
      - 2 \sqrt{\tfrac{5}{6}} \, \beta \sqrt{S_n m_\text{eff}} \, v_n,
\qquad S_n = 2 E^* \sqrt{r_\text{eff} \delta},

with overlap $\delta$, $v_n$ the relative velocity projected on the unit normal from $j$ to $i$,
the effective radius $r_\text{eff} = r_i r_j / (r_i + r_j)$ and mass
$m_\text{eff} = m_i m_j / (m_i + m_j)$ of the pair ($r_i$ and $m_i$ against a wall), the
effective modulus $E^* = E / (2 (1 - \nu^2))$ of two bodies of the same material (walls
included), and the damping ratio $\beta = -\ln e / \sqrt{\ln^2 e + \pi^2}$ for the coefficient of
restitution $e$, which the contact then reproduces nearly independently of the impact velocity.
The elastic energy of a contact is $\frac{8}{15} E^* \sqrt{r_\text{eff}} \, \delta^{5/2}$.

The parameters are `youngs_modulus` ($E$), `poissons_ratio` ($\nu$), and `restitution` ($e$, one
for no damping) of [KokkosParticleCloud.md]. The `stack` test checks the force law through the
equilibrium overlaps of two spheres stacked on a wall under gravity, and the `hertz_impact`
validation the force history, maximum overlap, and contact duration of an elastic impact against
Hertz's impact solution.

## Kuwabara-Kono damping

With a positive `dissipation_time` $A$, the damping is instead the viscoelastic one of Kuwabara
and Kono (1987) and Brilliantov et al. (1996), the dissipative constant times the rate of the
elastic force,

!equation
F_n = \frac{4}{3} E^* \sqrt{r_\text{eff}} \, \delta^{3/2} - \frac{3}{2} A \, \frac{4}{3} E^*
\sqrt{r_\text{eff} \delta} \, v_n,

whose coefficient of restitution is not a material constant but falls with the impact velocity,
$e = 1 - 1.15344 \, A \kappa^{2/5} v^{1/5} + \dots$ with $\kappa = (3/2)^{3/2} \cdot
\frac{4}{3} E^* \sqrt{r_\text{eff}} / m_\text{eff}$ (Schwager and Poeschel 1998; Ramirez et al.
1999). The `kk_*` validation tests check it against the integrated collision at four speeds
spanning a decade and a half. The total force may be attractive at the end of a contact, as in
the model; it is not clipped.
