#!/usr/bin/env python3
"""Generate the initial state and gold for periodic_gas.i (plan case V5): 64 spheres on a lattice
in a periodic unit box with pseudo-random velocities of zero total momentum and no damping, so the
total energy, linear momentum, and angular momentum are conserved. Prints the position and
velocity blocks to paste into the input and writes the gold of the conserved quantities."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

n_side, r, rho, v_max = 4, 0.03, 1000.0, 2.0
dt, num_steps = 0.01, 30
m = rho * 4.0 / 3.0 * math.pi * r**3

# Lattice positions and a fixed linear congruential sequence for the velocities, so the input is
# reproducible without a random module version dependence
positions, velocities = [], []
state = 12345
def lcg():
    global state
    state = (1103515245 * state + 12345) % 2**31
    return state / 2**31
for i in range(n_side):
    for j in range(n_side):
        for k in range(n_side):
            positions.append([(i + 0.5) / n_side, (j + 0.5) / n_side, (k + 0.5) / n_side])
            velocities.append([v_max * (2 * lcg() - 1) for _ in range(3)])
n = len(positions)
mean = [sum(v[c] for v in velocities) / n for c in range(3)]
for v in velocities:
    for c in range(3):
        v[c] = round(v[c] - mean[c], 6)
# Rounding leaves a residual momentum; take it out of the last particle exactly
for c in range(3):
    velocities[-1][c] -= sum(v[c] for v in velocities)

ke = sum(0.5 * m * sum(vc**2 for vc in v) for v in velocities)
print(f"m = {m:.6e}, N = {n}, KE = {ke:.12e}, volume fraction = {n * 4 / 3 * math.pi * r**3:.4f}")
print("initial_positions = '" + "\n                     ".join(" ".join(f"{x:g}" for x in p) for p in positions) + "'")
print("initial_velocities = '" + "\n                      ".join(" ".join(repr(x) for x in v) for v in velocities) + "'")

# Columns in the order MOOSE writes them; the ignored ones are fluctuating diagnostics
with open(os.path.join(here, "gold", "periodic_gas_out.csv"), "w") as f:
    f.write("time,angular_momentum,coordination_number,linear_momentum,load_imbalance,"
            "num_contacts,num_particles,total_energy,virial_pressure\n")
    for step in range(num_steps + 1):
        f.write(f"{step * dt:g},0,0,0,0,0,{n},{ke:.17g},0\n")
