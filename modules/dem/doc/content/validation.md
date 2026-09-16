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
