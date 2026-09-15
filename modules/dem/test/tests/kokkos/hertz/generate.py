#!/usr/bin/env python3
"""Generate the gold for stack.i (plan case V4): two equal spheres stacked on a wall under gravity
settle to the overlaps at which the Hertz force balances the weight above each contact,
  F = (4/3) E* sqrt(r_eff) delta^(3/2),
with r_eff = r against the wall (F = 2 m g) and r / 2 between the spheres (F = m g)."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

r, rho, E, nu, g = 0.01, 2500.0, 1e7, 0.3, 9.81
x = 0.05          # both spheres on this vertical line; the wall is y = 0
num_steps = 10
m = rho * 4.0 / 3.0 * math.pi * r**3
E_star = E / (2 * (1 - nu**2))

def overlap(force, r_eff):
    return (force / (4.0 / 3.0 * E_star * math.sqrt(r_eff)))**(2.0 / 3.0)

d_wall = overlap(2 * m * g, r)
d_pair = overlap(m * g, r / 2)
y0 = r - d_wall
y1 = y0 + 2 * r - d_pair
print(f"m = {m:.6e}, E* = {E_star:.6e}, wall overlap = {d_wall:.6e}, pair overlap = {d_pair:.6e}")
# Contact stiffness and period at equilibrium, to size the substep
for name, d, r_eff, m_eff in (("wall", d_wall, r, m), ("pair", d_pair, r / 2, m / 2)):
    k = 2 * E_star * math.sqrt(r_eff * d)
    print(f"{name}: stiffness {k:.4e} N/m, period {2 * math.pi * math.sqrt(m_eff / k):.4e} s")

with open(os.path.join(here, "gold", f"stack_out_state_{num_steps:04d}.csv"), "w") as f:
    f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
    f.write(f"0,0,0,0,0,0,0,{x},{y0:.17g},0\n")
    f.write(f"1,0,0,0,0,0,0,{x},{y1:.17g},0\n")
