#!/usr/bin/env python3
"""Generate the gold for corner.i (plan case V12, wedge): a sphere hits the reentrant corner of an
L-shaped mesh head-on, where two wall segments meet in a vertex convex toward the sphere. The
two segments report the same closest point, reduced to one contact along the line from the
vertex to the center, so the collision is the damped oscillator of a sphere against a wall
(effective mass m) with restitution e = exp(-zeta pi / sqrt(1 - zeta^2)): the velocity is
reversed with e and the center ends r + e v0 (T - t_touch - t_contact) from the vertex."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

r, rho, k, gamma = 0.05, 1000.0, 1e5, 64.7
x0, vertex, v0 = 0.8, 1.0, 1.0
dt, num_steps = 0.05, 8
T = dt * num_steps
m = rho * 4.0 / 3.0 * math.pi * r**3
omega0 = math.sqrt(k / m)
zeta = gamma / (2 * math.sqrt(k * m))
e = math.exp(-zeta * math.pi / math.sqrt(1 - zeta**2))
t_contact = math.pi / (omega0 * math.sqrt(1 - zeta**2))
gap = (vertex - x0) * math.sqrt(2) - r
t_touch = gap / v0
assert T > t_touch + t_contact
distance = r + e * v0 * (T - t_touch - t_contact)
c = 1 / math.sqrt(2)
print(f"e = {e:.9f}, contact from t = {t_touch:.6f} for {t_contact:.6e} s")
with open(os.path.join(here, "gold", f"corner_out_state_{num_steps:04d}.csv"), "w") as f:
    f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
    f.write(f"0,{-e * v0 * c:.17g},{-e * v0 * c:.17g},0,0,0,0,{vertex - distance * c:.17g},{vertex - distance * c:.17g},0\n")
