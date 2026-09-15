# LinearSpringDashpot

Linear spring-dashpot normal contact model used by [KokkosParticleCloud.md] (plan layer L4). It is
a compile-time policy struct with a `KOKKOS_INLINE_FUNCTION` force so the contact kernel inlines it
(plan decision D5). The normal force on particle $i$ from particle $j$ is

!equation
F_n = k_n \delta - \gamma_n v_n,

with overlap $\delta = r_i + r_j - |x_i - x_j|$ and $v_n$ the relative velocity projected on the
unit normal from $j$ to $i$. For two spheres of effective mass $m_\text{eff}$ the relative motion
is a damped linear oscillator, so the coefficient of restitution is
$e = \exp(-\zeta \pi / \sqrt{1 - \zeta^2})$ with $\zeta = \gamma_n / (2 \sqrt{k_n m_\text{eff}})$
and the contact lasts $\pi / (\omega_0 \sqrt{1 - \zeta^2})$, which the `restitution` test checks.

The parameters are `normal_stiffness` ($k_n$) and `normal_damping` ($\gamma_n$) of
[KokkosParticleCloud.md]; a zero stiffness disables contact. It is the default `contact_model`;
[Hertz.md] is the alternative.
