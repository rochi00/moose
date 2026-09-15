#!/usr/bin/env python3
"""Generate the analytic gold files for the integration tests (plan Section 8: references are
reproducible, not stored numbers)."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
gold = os.path.join(here, "gold")
os.makedirs(gold, exist_ok=True)

# --- ballistic.i: x = 0.1 + 2t, y = 0.9 + t - 5t^2, vx = 2, vy = 1 - 10t, dt = 0.04, 5 steps
for step in range(1, 6):
    t = 0.04 * step
    with open(os.path.join(gold, f"ballistic_out_state_{step:04d}.csv"), "w") as f:
        f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
        f.write(f"0,2,{1 - 10 * t:.17g},0,0,0,0,{0.1 + 2 * t:.17g},{0.9 + t - 5 * t * t:.17g},0\n")

# --- free_rotation.i: r = 0.1, rho = 1000, w = (1, 2, 3), 2 steps of dt = 1
r, rho = 0.1, 1000.0
m = rho * 4.0 / 3.0 * math.pi * r**3
I = 0.4 * m * r * r
w2 = 1 + 4 + 9
with open(os.path.join(gold, "free_rotation_out.csv"), "w") as f:
    f.write("time,angular_momentum,rotational_kinetic_energy\n")
    for t in range(3):
        f.write(f"{t},{I * math.sqrt(w2):.17g},{0.5 * I * w2:.17g}\n")
