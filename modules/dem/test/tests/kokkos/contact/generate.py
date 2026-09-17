#!/usr/bin/env python3
"""Generate the gold for restitution.i (plan case V3): two equal spheres in a head-on linear
spring-dashpot collision form a damped linear oscillator in their relative coordinate, so the
coefficient of restitution, contact duration, and final positions are analytic. Also the gold for
wall_restitution.i: the same oscillator for one sphere against a fixed wall, with the sphere's own
mass as the effective mass."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

r, rho = 0.05, 1000.0
m = rho * 4.0 / 3.0 * math.pi * r**3
k, gamma = 1e5, 64.7
v0 = 1.0          # each sphere approaches the origin at v0
x0 = 0.06         # spheres start at -x0 and +x0, so the gap is 2 (x0 - r)
dt, num_steps = 0.005, 4
T = dt * num_steps


def oscillator(m_eff):
    """Restitution and contact duration of the damped oscillator with the given effective mass"""
    omega0 = math.sqrt(k / m_eff)
    zeta = gamma / (2 * math.sqrt(k * m_eff))
    omega_d = omega0 * math.sqrt(1 - zeta**2)
    e = math.exp(-zeta * math.pi / math.sqrt(1 - zeta**2))
    t_touch = (x0 - r) / v0
    t_contact = math.pi / omega_d
    assert T > t_touch + t_contact
    print(f"m_eff = {m_eff:.6f}: restitution e = {e:.9f}, contact duration = {t_contact:.6e} s, "
          f"contact from t = {t_touch} to {t_touch + t_contact:.6f}")
    return e, t_touch, t_contact


e, t_touch, t_contact = oscillator(m / 2)

# The center of mass is at rest at the origin, so the relative coordinate carries everything:
# before contact the spheres approach at v0, after it they separate at e v0. The same solution
# holds along any direction; restitution.i collides along x in a 2D mesh and restitution_3d.i
# along the (1, 1, 1) diagonal of a 3D mesh, exercising every component.
xA = -r - e * v0 * (T - t_touch - t_contact)
# restitution_periodic.i collides the spheres the other way, through the periodic face at
# x = -0.2 and 0.2: sphere 0 ends at -0.2 - xA moving at +e v0 and sphere 1 mirrored.
# Each sphere ends at its contact point plus xA along its approach direction
diagonal = tuple(1 / math.sqrt(3) for _ in range(3))
for name, direction, contact_points in (
        ("restitution", (1, 0, 0), ((0, 0.5, 0), (0, 0.5, 0))),
        ("restitution_3d", diagonal, ((0, 0, 0), (0, 0, 0))),
        ("restitution_periodic", (-1, 0, 0), ((-0.2, 0.5, 0), (0.2, 0.5, 0)))):
    with open(os.path.join(here, "gold", f"{name}_out_state_{num_steps:04d}.csv"), "w") as f:
        f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
        for gid, sign in ((0, 1), (1, -1)):
            v = [f"{-sign * e * v0 * d:.17g}" for d in direction]
            x = [f"{contact_points[gid][c] + sign * xA * direction[c]:.17g}" for c in range(3)]
            f.write(f"{gid},{','.join(v)},0,0,0,{','.join(x)}\n")

# Against a fixed wall at x = 0 the sphere's own mass is the effective mass. wall_restitution.i
# hits the wall obliquely: the normal velocity component is reflected with e and the tangential
# one is unchanged since the contact is frictionless.
e, t_touch, t_contact = oscillator(m)
xA = -r - e * v0 * (T - t_touch - t_contact)
vt, y0 = 0.5, 0.3
with open(os.path.join(here, "gold", f"wall_restitution_out_state_{num_steps:04d}.csv"), "w") as f:
    f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
    f.write(f"0,{-e * v0:.17g},{vt:.17g},0,0,0,0,{xA:.17g},{y0 + vt * T:.17g},0\n")

# restitution_limit.i: the same collision with limit_damping. The relative coordinate follows
# the damped oscillator only until the dashpot's pull exceeds the spring's push on the way out,
# k delta = gamma |delta_dot| (delta the overlap, delta_dot < 0); from there the force is clamped
# at zero and the spheres coast apart at that velocity, which they keep after separating. So the
# restitution is higher than the free oscillator's, and the separation comes at t_touch + t* +
# delta(t*) / |delta_dot(t*)| instead of t_touch + pi / omega_d.
m_eff = m / 2
omega0 = math.sqrt(k / m_eff)
zeta = gamma / (2 * math.sqrt(k * m_eff))
omega_d = omega0 * math.sqrt(1 - zeta**2)
v_rel = 2 * v0                                  # closing speed of the relative coordinate
delta = lambda t: v_rel / omega_d * math.exp(-zeta * omega0 * t) * math.sin(omega_d * t)
delta_dot = lambda t: v_rel * math.exp(-zeta * omega0 * t) * (math.cos(omega_d * t) - zeta * omega0 / omega_d * math.sin(omega_d * t))
# Bisect k delta + gamma delta_dot = 0 on the rebound, between the turning point and separation
lo, hi = math.pi / (2 * omega_d), math.pi / omega_d
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if k * delta(mid) + gamma * delta_dot(mid) > 0:
        lo = mid
    else:
        hi = mid
t_star = 0.5 * (lo + hi)
e_limit = -delta_dot(t_star) / v_rel
t_touch = (x0 - r) / v0
t_end = t_touch + t_star + delta(t_star) / -delta_dot(t_star)
print(f"limit_damping: force clamped from t* = {t_star:.6e} s after touch, restitution {e_limit:.9f} "
      f"(free oscillator {math.exp(-zeta * math.pi / math.sqrt(1 - zeta**2)):.9f}), separation at t = {t_end:.6f}")
xA = -r - e_limit * v0 * (T - t_end)
with open(os.path.join(here, "gold", f"restitution_limit_out_state_{num_steps:04d}.csv"), "w") as f:
    f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
    for gid, sign in ((0, 1), (1, -1)):
        f.write(f"{gid},{-sign * e_limit * v0:.17g},0,0,0,0,0,{sign * xA:.17g},0.5,0\n")

# restitution_disk.i: the same collision with disk masses (density times pi r^2, LAMMPS's 2D
# convention): the damping ratio and contact time follow the lighter effective mass
m_disk = rho * math.pi * r**2
omega0 = math.sqrt(k / (m_disk / 2))
zeta = gamma / (2 * math.sqrt(k * m_disk / 2))
omega_d = omega0 * math.sqrt(1 - zeta**2)
e_disk = math.exp(-zeta * math.pi / math.sqrt(1 - zeta**2))
t_touch = (x0 - r) / v0
t_contact = math.pi / omega_d
print(f"disk masses: m = {m_disk:.6f}, restitution e = {e_disk:.9f}, contact duration = {t_contact:.6e} s")
# The heavier disks (7.85 kg against the spheres' 0.52) touch for 19.7 ms, past the 20 ms of
# the other cases, so this one runs 8 steps
T_disk = 8 * dt
assert T_disk > t_touch + t_contact
xA = -r - e_disk * v0 * (T_disk - t_touch - t_contact)
with open(os.path.join(here, "gold", "restitution_disk_out_state_0008.csv"), "w") as f:
    f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
    for gid, sign in ((0, 1), (1, -1)):
        f.write(f"{gid},{-sign * e_disk * v0:.17g},0,0,0,0,0,{sign * xA:.17g},0.5,0\n")
