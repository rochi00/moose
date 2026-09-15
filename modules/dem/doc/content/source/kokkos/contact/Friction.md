# Friction

Coulomb friction on the tangential spring-dashpot of the contact models of
[KokkosParticleCloud.md] (plan layer L4). The relative velocity of the contact point on particle
$i$ with respect to the partner, at $-r_i n$ from the center of $i$ and $r_j n$ from the center of
$j$ (the wall, rigid and at rest, has none),

!equation
v_\text{rel} = v_i - v_j - (r_i \omega_i + r_j \omega_j) \times n,

is split into its normal and tangential parts. The tangential spring displacement $\delta_t$ of
the contact, kept in its history ([KokkosParticleCloud.md#pair-states]), is projected on the
current tangent plane and stretched by $v_t \, \Delta t$ every substep, and the tangential force
on $i$ is

!equation
F_t = -k_t \delta_t - \gamma_t v_t, \qquad |F_t| \le \mu \max(F_n, 0),

with $k_t$ and $\gamma_t$ from the model ([LinearSpringDashpot.md]: `tangential_stiffness` and
`tangential_damping`; [Hertz.md]: the Mindlin no-slip stiffness $S_t = 8 G^* \sqrt{r_\text{eff}
\delta}$ and the Tsuji damping on it). When the limit is reached, the contact slides: the force
is scaled to the limit and the spring reset so that spring and dashpot together give exactly it,
as in LAMMPS's granular pair styles. The force acts at the contact point, so it exerts the torque
$-r_i n \times F_t$ on $i$; the partner gets the opposite force and, the arm being $r_j n$, the
torque $-r_j n \times F_t$.

The coefficient of friction is the `friction` parameter of [KokkosParticleCloud.md]; zero, or a
zero tangential stiffness, gives frictionless contact. The `sliding_sphere` test checks the
transition from slipping to rolling of a sphere launched along a wall against the analytic
$t^* = 2 v_0 / (7 \mu g)$ and rolling speed $5 v_0 / 7$.

## Rolling resistance

Rolling resistance is the elastic-plastic spring-dashpot model (Luding 2008; Ai et al. 2011,
model C; LAMMPS's `rolling sds`) on the rolling displacement $\delta_r$, the integral of the
rolling velocity $u = r_\text{eff} (\omega_i - \omega_j) \times n$ kept in the contact history
and projected on the tangent plane every substep. The resistance

!equation
F_r = -k_r \delta_r - \gamma_r u, \qquad |F_r| \le \mu_r \max(F_n, 0),

with the same reset on reaching the limit as the friction, acts as the torque
$r_\text{eff} \, n \times F_r$ on $i$ and its opposite on the partner, so the torque is limited
to $\mu_r F_n r_\text{eff}$. The rolling velocity, and with it $\delta_r$, is the same seen from
either side of the pair. The parameters are `rolling_friction` ($\mu_r$), `rolling_stiffness`
($k_r$), and `rolling_damping` ($\gamma_r$) of [KokkosParticleCloud.md]; zero rolling friction or
stiffness gives no rolling resistance. The `incline_rest`, `incline_roll`, and `incline_slide`
tests check a sphere on an inclined plane at rest, rolling under the plastic torque, and sliding
past the critical angle $\arctan \mu$ with its rotation locked, and `rolling_stop` the plastic
deceleration $\mu_r g / (1 + 2/5)$ of a rolling sphere.
