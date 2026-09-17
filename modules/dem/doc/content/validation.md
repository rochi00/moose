# DEM module validation

Beyond the verification tests of each object, the `kokkos/validation` tests reproduce the
contact-level results the DEM literature relies on. Each case is an input in
`modules/dem/test/tests/kokkos/validation` with a Python script generating its reference.

## Normal elastic impact (Hertz)

`hertz_impact`: two 10 mm spheres (E = 100 MPa, nu = 0.3, 2500 kg/m3) collide head-on at 1 m/s
with the elastic Hertz model. The force, position, and velocity histories, output every 10
microseconds over the 0.83 ms contact, match the integration of Hertz's impact equation
$m_\text{eff} \ddot\delta = -\frac{4}{3} E^* \sqrt{r_\text{eff}} \delta^{3/2}$ to 1.2e-4; the
maximum overlap is within 1.4e-4 of $(15 m_\text{eff} v^2 / (16 E^* \sqrt{r_\text{eff}}))^{2/5}$
and the contact lasts 0.81-0.83 ms against $2.94 (m_\text{eff}^2 / (r_\text{eff} E^{*2} v))^{1/5}
= 0.832$ ms.

## Coefficient of restitution vs impact velocity

`kk_v01` to `kk_v3`: the same spheres with the Kuwabara-Kono viscoelastic damping
(`dissipation_time = 5e-5` s) at closing speeds of 0.1, 0.3, 1, and 3 m/s. The rebound matches
the integrated viscoelastic collision to 1e-4 (3e-4 at 3 m/s); e falls from 0.809 to 0.661, with
$1 - e$ following the $v^{1/5}$ law at low speed (the first-order series of Schwager and Poeschel
is within 1-4%), where the Tsuji damping keeps e at its set value (the `restitution` tests).
This is the comparison of Thornton, Cummins, and Cleary (2013) between constant and
velocity-dependent restitution models.

## Oblique impact

`oblique_60` and `oblique_impact_sweep.py`: a 5 mm aluminium oxide sphere impacting a soda-lime
glass plane at 3.85 m/s with $e_n = 0.98$ and $\mu = 0.092$, the setting of Kharaz, Gorham, and
Salman (2001); the cloud's single material is chosen to give the pair's effective moduli
($E^* = 63$ GPa, $G^* = 13.5$ GPa). Above the critical angle
$\arctan(\tfrac{7}{2} \mu (1 + e_n)) = 32.5^\circ$ the contact slides throughout and the rebound
spin $\tfrac{5}{2} \mu (1 + e_n) v_n / r$ and tangential velocity $v_t - \mu (1 + e_n) v_n$ are
analytic, matched by the test to 2e-3 at $60^\circ$ and by the sweep at every angle above it.
Below it the contact sticks for part of the impact and the rebound depends on the tangential
model, as in the theory of Maw, Barber, and Fawcett (1976): the sweep shows the rebound spin
peaking at about 610 rad/s just below the critical angle and the tangential restitution
$-(v_t' + r\omega') / v_t$ turning negative at small angles (the contact point springs back),
with the no-slip and partial-slip springs differing only there.

!media media/dem_oblique_impact_sweep.png
       caption=Rebound spin and tangential restitution against the impact angle for the two
               tangential springs, with the sliding-throughout solution.

## Hertzian chain: Nesterenko's solitary wave

`nesterenko` and `nesterenko_check.py`: sixty 10 mm steel beads (Y = 203 GPa, nu = 0.3,
7780 kg/m3) barely touching in a line with the elastic Hertz contact; the first bead is set
moving at 1 m/s and a single solitary wave forms and travels the chain (Nesterenko 1983; Job,
Melo, Sokolow, and Sen 2008, Eq. 2). Measured from the bead velocities output every 2
microseconds:

- the wave speed at the amplitude the wave settles to (0.68 m/s) is within 0.8% of the long
  wavelength solution $v = (6 / (5 \pi \rho \theta))^{2/5} V_m^{1/5}$, $\theta = 3 (1 - \nu^2) /
  (4 Y)$, and is constant to 1e-5 along the chain, as is the amplitude: a genuine solitary wave;
- the pulse spans 3.96 bead diameters at 1% of its peak against 3.95 for the $\cos^4$ profile of
  characteristic length $R \sqrt{10}$ (whose full support is 5 diameters), while its half-maximum
  width is 13% narrower than $\cos^4$: the continuum profile's neglected terms are of order
  $(2R/\lambda)^5 = 10\%$, the discrete wave being sharper in its core (Chatterjee 1999).

## Homogeneous cooling: Haff's law

`haff` and `haff_check.py`: 512 frictionless spheres at a solid fraction of 2% in a periodic
cube, linear spring-dashpot contact set for $e = 0.9$, Gaussian initial velocities, no walls, no
gravity. Kinetic theory (Haff 1983; Brilliantov and Poeschel 2004) gives
$T(t) = T_0 / (1 + t/\tau_0)^2$ with $1/\tau_0 = (1 - e^2)\,\nu(T_0)/6$ for the Enskog collision
frequency $\nu = 4\sqrt{\pi} n \sigma^2 g_2(\phi) \sqrt{T/m}$. Over four cooling times the fitted
$\tau_0$ is within 3% of the prediction (0.752 s against 0.731 s) and $1/\sqrt{T/T_0}$ stays
within 2% of a straight line; $T$ falls to 0.040 $T_0$ against Haff's 0.038. Setting this case up
exposed and fixed a limitation of the periodic and partition handling: a particle straying more
than the neighbor-list skin beyond its rank's box within one MOOSE step lost its far-side
neighbors, so such particles are now settled (wrapped, located, migrated) within the substep loop.

## Cross-code check against LAMMPS

`kokkos/lammps` (with `generate.py` and the LAMMPS inputs and dumps under `lammps/`): the
32-sphere configuration of LAMMPS's own granular unit tests (a jittered fcc lattice of spheres
of diameter 1.25 in a fully periodic box of side 3.36, every sphere overlapping its neighbours)
run for 100 steps in LAMMPS (22 Jul 2025) with `pair gran/hooke/history 1000 500 50 50 mu 1`
and in this module with the `LinearSpringDashpot` model, the LAMMPS damping coefficients times
$m_\text{eff}$; and a spinning sphere striking a `fix wall/gran` floor obliquely against an
analytic wall. Both codes integrate with velocity Verlet and the same contact conventions, and
the LAMMPS dumps are the gold of the tests: positions agree to $5 \times 10^{-14}$ and velocities
and spins to $10^{-14}$ over the 100 steps, frictionless and with the contacts sliding at the
Coulomb limit (the LAMMPS-style spring reset), and through the loaded phase of the wall bounce.
Two conventions differ and are kept out of the reference windows: LAMMPS caps the tangential
force at $\mu |F_n|$ where this module uses $\mu \max(F_n, 0)$ (no friction without contact
pressure), which matters only when the dashpot makes the net normal force attractive; and
LAMMPS's setup-step force evaluation does not advance the histories, so a contact whose
dashpot alone exceeds the Coulomb limit at step 0 gets no tangential force on that step. A
`limit_damping` option (LAMMPS's, clamping the normal force at zero) is provided for
comparisons with `pair granular`. The same collision as `examples/granular/in.restitution`
(Hertz, `damping coeff_restitution`) rebounds at 0.79999999 m/s in both codes.

The comparison found and fixed one defect: the forces used for the first half-kick of a MOOSE
step were recomputed at the end of the previous step with the full-step velocities, while every
substep, like LAMMPS, uses forces computed with the half-step velocities, so the trajectory
depended on how the substeps were grouped into steps at the order of the dashpot force times
the substep. The last substep's forces now travel with the particles (and into checkpoints),
and the end-of-step recompute serves the outputs only (`friction/oblique_impact_regroup`).

## Dense-flow benchmarks (`examples/`)

The Tier 2 cases below are run by scripts under `modules/dem/examples/` that generate the
inputs, run them (8 MPI ranks on a laptop, 30 s to 6 min per case), and analyze the CSV output.
They are not regression tests: each is a few 10^4 spheres with linear spring-dashpot contacts
(overlaps below 1% of the radius, restitution 0.45, substep a 26th of the contact time) and the
result depends on the friction coefficients, which are the calibration parameters of DEM.

### Silo discharge: Beverloo's law

`silo_discharge_3d/run.py`: a cylindrical silo 30 $d$ wide filled to 40 $d$ with 30,314 beads
($d$ = 6 mm, radii spread by 5%) drains through a circular floor orifice of diameter $D$ = 6, 7, 8
and 10 $d$ (the inner circle of a `ConcentricCircleMeshGenerator` disk extruded along $z$, so the
orifice area is exact); beads through the orifice leave the mesh and are counted by `num_exited`.
Beverloo's law $W = C \rho_b \sqrt{g} (D - k d)^{5/2}$ with $C = 0.58$ and $k = 1.4$ (Beverloo et
al. 1961; Nedderman et al. 1982) is tested by fitting $W^{2/5}$ against $D$, a straight line whose
intercept is $k d$:

| friction | $C$ | $k$ | exponent |
|---|---|---|---|
| $\mu = 0.5$, $\mu_r = 0.1$ | 0.49 | 1.68 | 2.47 |
| $\mu = 0.3$, $\mu_r = 0$ (glass-bead-like) | 0.548 | 1.57 | 2.49 |

Every rate is steady to better than 1% over its window (free surface above one silo diameter) and
independent of the fill height; the measured packing fraction is 0.606. The functional form is
exact; the prefactor is within 6% of Beverloo's with glass-bead-like friction and 15% low with
rolling resistance, the direction and size reported for rolling resistance in the literature.

!media media/dem_silo_discharge_3d.png style=width:60% caption=Discharge rate against orifice diameter, glass-bead-like friction.

### Axisymmetric column collapse

`column_collapse_3d/run.py`: 22,801 beads rained into a cylinder of radius $R_0 = 12 d$ and
settled (packing 0.58), sliced at $H_0 = a R_0$ and released on a flat floor for $a$ = 0.5, 1, 2, 3
($\mu = 0.5$, $\mu_r = 0.1$). The deposit edge (outermost annulus of width $d$ still covered by
30% of a monolayer) against Lube et al. 2004, $(R_\infty - R_0)/R_0 = 1.24 a$ for $a < 1.7$ and
$1.6 a^{1/2}$ above, and Lajeunesse et al. 2004, $1.0 a$ / $2.0 a^{1/2}$ about $a = 0.74$:

| $a$ | DEM | Lube | Lajeunesse |
|---|---|---|---|
| 0.5 | 0.50 | 0.62 | 0.50 |
| 1 | 1.25 | 1.24 | 2.0 |
| 2 | 2.42 | 2.26 | 2.83 |
| 3 | 3.42 | 2.77 | 3.46 |

The runouts lie within the band spanned by the two experiments (rough grains for Lube, glass
beads for Lajeunesse), with the deposits at rest (kinetic energy at $10^{-4}$ of its peak by
$2 t_\infty$). The deposit heights for $a \ge 1$, 0.63 to 0.71 $R_0$, are below Lube's ~1 $R_0$:
smooth spheres with $\mu_r = 0.1$ leave a shallow cone where rough grains keep a steep core.

!media media/dem_column_collapse_3d.png style=width:80% caption=Runout and deposit height against aspect ratio.

### Angle of repose

`angle_of_repose/run.py` reproduces the setup of Zhou, Xu, Yu and Zulli 2002 (Powder Technol.
125): 10 mm spheres settled on a 30 $d$ plate in a container 40 $d$ wide and $w$ thick are
discharged over both edges of the plate; the angle is the slope of the heap's surface profile
between 15% and 85% of the plate half-width, averaged over both flanks. Sliding friction as
their base condition ($\mu_{pp} = 0.4$, $\mu_{pw} = 0.6$, using `wall_friction`). Their rolling
model is a constant torque against each particle's own spin, so their rolling coefficient does
not carry over to the spring-dashpot model on relative rolling used here (mapped as if it did,
their base $\mu_r$ gives 18.7 degrees against their 28.5); the test is therefore against their
glass-bead experiments, $\theta_0 = 29.6 d^{-0.206}$ degrees without front and rear walls (18.4
for 10 mm) times $1 + e^{-0.18 w/d}$ with them (27.4 at $w = 4 d$, 22.8 at $8 d$):

| $\mu_r$ | $w = 4 d$ | $w = 8 d$ | periodic ($w \to \infty$) |
|---|---|---|---|
| 0.02 | 18.7 | | |
| 0.05 | 29.2 | 27.2 | 19.4 |
| 0.1 | 37.2 | 33.3 | 27.2 |
| 0.2 | 47.0 | | |
| experiment | 27.4 | 22.8 | 18.4 |

With a single rolling coefficient of 0.05 the heap angle matches the experiment within 2 degrees
(the paper's own measurement error) at $w = 4 d$ and without walls, and the wall-effect ratio
$\theta(4d)/\theta_0$ is 1.51 against 1.49 measured; the $8 d$ case is 4 degrees high. Every heap
comes to rest (kinetic energy below $10^{-14}$ J), and the angle rises monotonically with
$\mu_r$, both flanks agreeing to 1-4 degrees on heaps of 300-1200 spheres.

!media media/dem_angle_of_repose.png style=width:80% caption=Angle of repose against rolling friction coefficient at w = 4 d, and against container thickness.
