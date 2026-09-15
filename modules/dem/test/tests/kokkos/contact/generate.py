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
# before contact the spheres approach at v0, after it they separate at e v0. The same solution
# holds along any direction; restitution.i collides along x in a 2D mesh and restitution_3d.i
# along the (1, 1, 1) diagonal of a 3D mesh, exercising every component.
T = dt * num_steps
assert T > t_touch + t_contact
xA = -r - e * v0 * (T - t_touch - t_contact)
for name, direction, center in (("restitution", (1, 0, 0), (0, 0.5, 0)),
                                ("restitution_3d", tuple(1 / math.sqrt(3) for _ in range(3)), (0, 0, 0))):
    with open(os.path.join(here, "gold", f"{name}_out_state_{num_steps:04d}.csv"), "w") as f:
        f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
        for gid, sign in ((0, 1), (1, -1)):
            v = [f"{-sign * e * v0 * d:.17g}" for d in direction]
            x = [f"{center[c] + sign * xA * direction[c]:.17g}" for c in range(3)]
            f.write(f"{gid},{','.join(v)},0,0,0,{','.join(x)}\n")
