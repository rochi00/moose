#!/usr/bin/env python3
"""Fit Haff's law to the kinetic energy history of haff.i and compare the cooling time with
kinetic theory: T(t) = T0 / (1 + t / tau0)^2, 1 / tau0 = (1 - e^2) nu(T0) / 6 (see
generate_haff.py). A least-squares line through 1 / sqrt(T / T0) against t gives the measured
1 / tau0; the fit also reports how far the history departs from the law. Exits nonzero when the
cooling time is off by more than the tolerance. Usage: haff_check.py [directory]"""
import csv, math, os, sys

here = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
L, m, T0_gen, nu0_gen, tau0_gen, e, gamma = [float(x) for x in open(os.path.join(here, "haff_params.txt")).read().split()]
N, r, phi = 512, 0.005, 0.02

rows = list(csv.DictReader(open(os.path.join(here, "haff_out.csv"))))
t = [float(q["time"]) for q in rows]
ke = [float(q["kinetic_energy"]) for q in rows]
T = [2 * k / (3 * N) for k in ke]           # m <v^2> / 3 = (2/3) KE / N
T0 = T[0]
sigma, n = 2 * r, N / L**3
g2 = (2 - phi) / (2 * (1 - phi)**3)
nu0 = 4 * math.sqrt(math.pi) * n * sigma**2 * g2 * math.sqrt(T0 / m)
tau0 = 6 / ((1 - e**2) * nu0)

# Least squares of y = 1 / sqrt(T / T0) = 1 + t / tau0 through the origin offset 1
y = [1 / math.sqrt(Ti / T0) for Ti in T]
sxx = sum(ti * ti for ti in t)
sxy = sum(ti * (yi - 1) for ti, yi in zip(t, y))
slope = sxy / sxx
tau0_fit = 1 / slope
residual = max(abs(yi - (1 + ti * slope)) / yi for ti, yi in zip(t, y))
print(f"T0 = {T0 / m:.4f} m (m2/s2), Enskog collision frequency nu0 = {nu0:.2f} /s, "
      f"theory tau0 = {tau0:.4f} s; fit tau0 = {tau0_fit:.4f} s, deviation {100 * abs(tau0_fit - tau0) / tau0:.2f}%, "
      f"largest departure of 1/sqrt(T/T0) from the line {100 * residual:.2f}%")
print(f"T fell to {T[-1] / T0:.4f} T0 at t = {t[-1]} s (Haff: {1 / (1 + t[-1] / tau0)**2:.4f})")
if abs(tau0_fit - tau0) / tau0 > 0.06 or residual > 0.03:
    print("FAILED")
    sys.exit(1)
print("within tolerance")
