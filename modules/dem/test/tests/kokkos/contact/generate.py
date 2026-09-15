#!/usr/bin/env python3
"""Generate the gold for restitution.i (plan case V3): two equal spheres in a head-on linear
spring-dashpot collision form a damped linear oscillator in their relative coordinate, so the
coefficient of restitution, contact duration, and final positions are analytic."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

r, rho = 0.05, 1000.0
m = rho * 4.0 / 3.0 * math.pi * r**3
m_eff = m / 2
k, gamma = 1e5, 64.7
v0 = 1.0          # each sphere approaches the origin at v0
x0 = 0.06         # spheres start at -x0 and +x0, so the gap is 2 (x0 - r)
dt, num_steps = 0.005, 4

omega0 = math.sqrt(k / m_eff)
zeta = gamma / (2 * math.sqrt(k * m_eff))
omega_d = omega0 * math.sqrt(1 - zeta**2)
e = math.exp(-zeta * math.pi / math.sqrt(1 - zeta**2))
t_touch = (x0 - r) / v0
t_contact = math.pi / omega_d
print(f"restitution e = {e:.9f}, contact duration = {t_contact:.6e} s, "
      f"contact from t = {t_touch} to {t_touch + t_contact:.6f}")

# The center of mass is at rest at the origin, so the relative coordinate carries everything:
# before contact the spheres approach at v0, after it they separate at e v0
with open(os.path.join(here, "gold", "restitution_out_state_%04d.csv" % num_steps), "w") as f:
    T = dt * num_steps
    assert T > t_touch + t_contact
    xA = -r - e * v0 * (T - t_touch - t_contact)
    f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
    f.write(f"0,{-e * v0:.17g},0,0,0,0,0,{xA:.17g},0.5,0\n")
    f.write(f"1,{e * v0:.17g},0,0,0,0,0,{-xA:.17g},0.5,0\n")
