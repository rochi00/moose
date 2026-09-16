#!/usr/bin/env python3
"""Generate the gold for hertz_impact.i (Tier 1: normal elastic Hertz impact): two equal spheres
of radius r collide head-on at closing speed 2 v0 with the elastic Hertz model. The relative
motion obeys m_eff d'' = -(4/3) E* sqrt(r_eff) d^(3/2) from touch until separation, so the
maximum overlap and the contact duration are (Hertz 1882; Johnson, Contact Mechanics, 11.4)
  d_max = (15 m_eff v^2 / (16 E* sqrt(r_eff)))^(2/5),   t_c = 2.94 (m_eff^2 / (r_eff E*^2 v))^(1/5)
with v = 2 v0 the closing speed, and the force history F(t) follows from the same equation,
integrated here with fourth-order Runge-Kutta at a step far below the tolerance. The input
outputs the state every 1e-5 s, about 80 points over the contact; the gold holds the force,
position, and velocity of both spheres at eight of them spanning the contact, and the full
reference history is written next to it for plotting."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

r, rho, E, nu = 0.01, 2500.0, 1e8, 0.3
v0, x0 = 0.5, 0.01001           # each sphere approaches at v0 from +-x0: the gap is 2 (x0 - r)
dt, num_steps = 1e-5, 90         # the state is output every step; the contact spans ~83 of them
compared = (4, 20, 40, 50, 60, 70, 80, 90)
m = rho * 4.0 / 3.0 * math.pi * r**3
m_eff, r_eff = m / 2, r / 2
E_star = E / (2 * (1 - nu**2))
k = 4.0 / 3.0 * E_star * math.sqrt(r_eff)
v = 2 * v0
t_touch = (x0 - r) / v0
d_max = (15 * m_eff * v**2 / (16 * E_star * math.sqrt(r_eff)))**0.4
t_c = 2.94 * (m_eff**2 / (r_eff * E_star**2 * v))**0.2
print(f"m = {m:.6e}, E* = {E_star:.6e}, touch at t = {t_touch:.6e}, d_max = {d_max:.6e}, "
      f"t_c = {t_c:.6e}, F_max = {k * d_max**1.5:.6f}; the run covers t = {dt * num_steps:.2e}")
assert dt * num_steps > t_touch + t_c

# RK4 of (d, u) with u = d', the overlap and its rate; before touch d < 0 is minus the gap,
# closing at v
h = 1e-8
def rhs(state):
    d, u = state
    return (u, -k * max(d, 0)**1.5 / m_eff)
state = (-2 * (x0 - r), v)
t = 0.0
with open(os.path.join(here, "gold", "hertz_impact_reference.csv"), "w") as ref:
    ref.write("time,force,overlap\n")
    for step in range(1, num_steps + 1):
        while t < step * dt - 0.5 * h:
            k1 = rhs(state)
            k2 = rhs(tuple(s + 0.5 * h * q for s, q in zip(state, k1)))
            k3 = rhs(tuple(s + 0.5 * h * q for s, q in zip(state, k2)))
            k4 = rhs(tuple(s + h * q for s, q in zip(state, k3)))
            state = tuple(s + h / 6 * (a + 2 * b + 2 * c + e) for s, a, b, c, e in zip(state, k1, k2, k3, k4))
            t += h
        d, u = state
        force = k * max(d, 0)**1.5
        ref.write(f"{step * dt:.17g},{force:.17g},{max(d, 0):.17g}\n")
        if step in compared:
            # Sphere 0 approaches from -x at v0 and is pushed back along -x; the center is at
            # -(r - d/2) from the contact plane, half the closing rate away from it
            with open(os.path.join(here, "gold", f"hertz_impact_out_state_{step:04d}.csv"), "w") as f:
                f.write("fx,fy,fz,gid,taux,tauy,tauz,vx,vy,vz,wx,wy,wz,x,y,z\n")
                f.write(f"{-force:.17g},0,0,0,0,0,0,{u / 2:.17g},0,0,0,0,0,{-(r - d / 2):.17g},0.5,0\n")
                f.write(f"{force:.17g},0,0,1,0,0,0,{-u / 2:.17g},0,0,0,0,0,{r - d / 2:.17g},0.5,0\n")
