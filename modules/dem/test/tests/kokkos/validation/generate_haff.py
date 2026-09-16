#!/usr/bin/env python3
"""Generate the initial state of haff.i (Tier 2: homogeneous cooling state) and print the
parameters of Haff's law for it. 512 inelastic frictionless spheres on a lattice in a periodic
cube at a solid fraction of 2%, with Gaussian random velocities of zero total momentum, cool by
their collisions alone. With a constant coefficient of restitution e (the linear spring-dashpot
model), kinetic theory gives (Haff 1983; Brilliantov and Poeschel 2004)
  dT/dt = -(1 - e^2) nu(T) T / 3,  nu(T) = 4 sqrt(pi) n sigma^2 g2(phi) sqrt(T / m),
for the granular temperature T = m <v^2> / 3 and the Enskog collision frequency nu with
g2 = (2 - phi) / (2 (1 - phi)^3), hence
  T(t) = T0 / (1 + t / tau0)^2,  1 / tau0 = (1 - e^2) nu(T0) / 6.
haff_check.py fits 1 / sqrt(T / T0) = 1 + t / tau0 to the kinetic energy history and compares."""
import math, os, random

here = os.path.dirname(os.path.abspath(__file__))
n_side, r, rho, phi, e, k = 8, 0.005, 2500.0, 0.02, 0.9, 1e5
N = n_side**3
m = rho * 4.0 / 3.0 * math.pi * r**3
L = (N * 4.0 / 3.0 * math.pi * r**3 / phi)**(1.0 / 3.0)
# Dashpot for the restitution e of a pair (effective mass m / 2)
zeta = -math.log(e) / math.sqrt(math.pi**2 + math.log(e)**2)
gamma = 2 * zeta * math.sqrt(k * m / 2)
sigma = 2 * r
n = N / L**3
g2 = (2 - phi) / (2 * (1 - phi)**3)

random.seed(7)
positions, velocities = [], []
for i in range(n_side):
    for j in range(n_side):
        for l in range(n_side):
            positions.append(((i + 0.5) * L / n_side, (j + 0.5) * L / n_side, (l + 0.5) * L / n_side))
            velocities.append([random.gauss(0, 1.5) for _ in range(3)])
mean = [sum(v[c] for v in velocities) / N for c in range(3)]
for v in velocities:
    for c in range(3):
        v[c] -= mean[c]
T0 = m * sum(sum(vc**2 for vc in v) for v in velocities) / (3 * N)
nu0 = 4 * math.sqrt(math.pi) * n * sigma**2 * g2 * math.sqrt(T0 / m)
tau0 = 6 / ((1 - e**2) * nu0)
t_c = math.pi / math.sqrt(k / (m / 2) * (1 - zeta**2))
print(f"N = {N}, L = {L:.6f}, m = {m:.4e}, normal_damping = {gamma:.6f}, contact {t_c:.2e} s, "
      f"T0/m = {T0 / m:.4f}, nu0 = {nu0:.3f} /s (mean free time {1 / nu0:.4f} s), tau0 = {tau0:.4f} s")

with open(os.path.join(here, "haff_init.txt"), "w") as f:
    f.write("initial_positions = '" + "\n                     ".join(" ".join(f"{x:.6g}" for x in p) for p in positions) + "'\n")
    f.write("initial_velocities = '" + "\n                      ".join(" ".join(repr(x) for x in v) for v in velocities) + "'\n")
with open(os.path.join(here, "haff_params.txt"), "w") as f:
    f.write(f"{L!r} {m!r} {T0!r} {nu0!r} {tau0!r} {e} {gamma!r}\n")
