#!/usr/bin/env python3
"""Generate the golds for kuwabara_kono.i (Tier 1: coefficient of restitution vs impact velocity):
two equal spheres collide head-on with the Hertz elastic force and the Kuwabara-Kono (1987) /
Brilliantov et al. (1996) viscoelastic damping, F = k (d^(3/2) + (3/2) A d^(1/2) d'), with
k = (4/3) E* sqrt(r_eff). The restitution coefficient is not a constant but falls with the
impact velocity: to first order (Schwager and Poeschel 1998; Ramirez et al. 1999)
  e = 1 - C1 A kappa^(2/5) v^(1/5) + ...,  kappa = (3/2)^(3/2) k / m_eff,  C1 = 1.15344,
which is printed for comparison. The reference used for the golds is the numerical integration
of the collision itself (fourth-order Runge-Kutta, until the overlap returns to zero), giving e
and the contact duration for each closing speed, from which the final velocities and positions
follow as in the restitution tests. Four closing speeds spanning a decade and a half check the
v^(1/5) trend."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

r, rho, E, nu, A = 0.01, 2500.0, 1e8, 0.3, 5e-5
m = rho * 4.0 / 3.0 * math.pi * r**3
m_eff, r_eff = m / 2, r / 2
E_star = E / (2 * (1 - nu**2))
k = 4.0 / 3.0 * E_star * math.sqrt(r_eff)
kappa = 1.5**1.5 * k / m_eff
gap = 2e-5                       # each sphere starts gap / 2 from touching: x0 = r + gap / 2
cases = {"kk_v01": 0.1, "kk_v03": 0.3, "kk_v1": 1.0, "kk_v3": 3.0}   # closing speeds v = 2 v0
dt, num_steps = 2e-4, 12
T = dt * num_steps


def collide(v):
    """Integrate one collision at closing speed v; return (e, t_contact)"""
    h = 2e-9
    d, u, t = 0.0, v, 0.0
    def acc(d, u):
        return -(k / m_eff) * (max(d, 0)**1.5 + 1.5 * A * math.sqrt(max(d, 0)) * u)
    while True:
        k1 = (u, acc(d, u))
        k2 = (u + 0.5 * h * k1[1], acc(d + 0.5 * h * k1[0], u + 0.5 * h * k1[1]))
        k3 = (u + 0.5 * h * k2[1], acc(d + 0.5 * h * k2[0], u + 0.5 * h * k2[1]))
        k4 = (u + h * k3[1], acc(d + h * k3[0], u + h * k3[1]))
        d_new = d + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        u_new = u + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        t += h
        if d_new <= 0 and t > h:
            # Linear interpolation of the exit time
            frac = d / (d - d_new)
            return -(u + frac * (u_new - u)) / v, t - h + frac * h
        d, u = d_new, u_new


for name, v in cases.items():
    e, t_c = collide(v)
    e1 = 1 - 1.15344 * A * kappa**0.4 * v**0.2
    v0 = v / 2
    t_touch = gap / v
    assert T > t_touch + t_c
    x = r + e * v0 * (T - t_touch - t_c)
    print(f"{name}: v = {v}: e = {e:.6f} (first-order series {e1:.6f}, 1 - e = {1 - e:.4e}), "
          f"contact {t_c:.4e} s, t_touch = {t_touch:.2e}")
    with open(os.path.join(here, "gold", f"{name}_state_{num_steps:04d}.csv"), "w") as f:
        f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
        f.write(f"0,{-e * v0:.17g},0,0,0,0,0,{-x:.17g},0.5,0\n")
        f.write(f"1,{e * v0:.17g},0,0,0,0,0,{x:.17g},0.5,0\n")
print(f"initial position +-{r + gap / 2!r}")
