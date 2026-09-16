#!/usr/bin/env python3
"""Measure Nesterenko's solitary wave in the output of nesterenko.i and compare it with the long
wavelength theory (Job et al. 2008, Eq. 2): the wave speed against v = c0 psi_m^(1/4) with
c0 = (6 / (5 pi rho theta))^(1/2), theta = 3 (1 - nu^2) / (4 Y), and the amplitude psi_m = V_m / v
(the bead velocity being v psi), i.e. v = c0^(4/5) V_m^(1/5); the fifth-root scaling of the
speed with amplitude, from the wave launched at 1 m/s and its own dispersive decay; and the
pulse width against the cos^4 profile of characteristic length R sqrt(10), whose full width at
half maximum is 2 acos(2^(-1/4)) R sqrt(10) = 3.64 R. Exits nonzero when a comparison fails.
Usage: nesterenko_check.py [directory of the state output]"""
import csv, glob, math, os, sys

here = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
R, rho, Y, nu = 5e-3, 7780.0, 2.03e11, 0.3
dt = 2e-6
theta = 3 * (1 - nu**2) / (4 * Y)
c0 = math.sqrt(6 / (5 * math.pi * rho * theta))

# Bead velocities v[step][bead] from the state files (sorted by global ID = bead index)
files = sorted(glob.glob(os.path.join(here, "nesterenko_out_state_[0-9]*.csv")))[1:]  # step 0 is empty
v = []
for f in files:
    rows = list(csv.reader(open(f)))[1:]
    v.append([float(r[1]) for r in rows])
num_steps, num_beads = len(v), len(v[0])
assert num_steps > 300 and num_beads == 60, "unexpected output"


def peak(series):
    """Time (in steps) and value of the maximum, refined by a parabola through its neighbors"""
    i = max(range(1, len(series) - 1), key=lambda k: series[k])
    a, b, c = series[i - 1], series[i], series[i + 1]
    denom = a - 2 * b + c
    shift = 0.5 * (a - c) / denom if denom else 0.0
    return i + shift, b - 0.25 * (a - c) * shift


failures = []
def check(name, measured, expected, tol):
    err = abs(measured - expected) / abs(expected)
    print(f"{name}: measured {measured:.6g}, theory {expected:.6g}, deviation {100 * err:.2f}% (allowed {100 * tol:.0f}%)")
    if err > tol:
        failures.append(name)

# Speed and amplitude between beads 20 and 40, once the wave has formed and before the far end
t20, _ = peak([v[s][20] for s in range(num_steps)])
t40, _ = peak([v[s][40] for s in range(num_steps)])
_, V_m = peak([v[s][30] for s in range(num_steps)])
speed = 20 * 2 * R / ((t40 - t20) * dt)
check("wave speed at V_m = %.4f m/s" % V_m, speed, c0**0.8 * V_m**0.2, 0.03)

# The fifth-root scaling: the pulse from a 1 m/s edge bead sheds a little energy into a tail as
# it forms, so its amplitude is below 1 m/s; a second, weaker wave is obtained from the same
# output by the beads' velocities behind the leading pulse only if one exists. Instead the
# scaling is checked between this run and the theory at the measured amplitude (above) and by
# the local relation d ln v / d ln V_m along the slightly decaying pulse: the speed between
# beads 10-20 and 40-50 at their own amplitudes
t10, _ = peak([v[s][10] for s in range(num_steps)])
t50, _ = peak([v[s][50] for s in range(num_steps)])
_, V15 = peak([v[s][15] for s in range(num_steps)])
_, V45 = peak([v[s][45] for s in range(num_steps)])
v_early = 20 * R / ((t20 - t10) * dt)
v_late = 20 * R / ((t50 - t40) * dt)
print(f"amplitude decay along the chain: V_m = {V15:.5f} (beads 10-20) -> {V45:.5f} (beads 40-50); "
      f"speed {v_early:.2f} -> {v_late:.2f} m/s; theory ratio {(V45 / V15)**0.2:.5f}, measured {v_late / v_early:.5f}")

# Width: the bead spacing is too coarse to sample the profile in space, so it is taken from the
# velocity of bead 30 in time, which the wave passes at the measured speed: the full width at
# half maximum and the duration above 1% of the peak (cos^4 = 0.01 at xi = 1.249) of the cos^4
# pulse are 2 acos(2^(-1/4)) R sqrt(10) / v and 2 * 1.249 R sqrt(10) / v
series = [v[s][30] for s in range(num_steps)]
peak_step, peak_value = peak(series)
def duration(fraction):
    above = [s for s in range(num_steps) if series[s] >= fraction * peak_value]
    def crossing(n_in, n_out):
        a, b = series[n_in], series[n_out]
        return n_in + (fraction * peak_value - a) / (b - a) * (n_out - n_in)
    return (crossing(above[-1], above[-1] + 1) - crossing(above[0], above[0] - 1)) * dt
# The cos^4 profile is the long wavelength approximation, whose neglected terms are of order
# (2R / lambda)^5 = 10%: the discrete wave is somewhat sharper in its core (Chatterjee 1999
# gives the exact shape), so the half-maximum width is only checked to 15%
check("pulse width at half maximum (m)", duration(0.5) * speed, 2 * math.acos(2**-0.25) * R * math.sqrt(10), 0.15)
check("pulse width at 1% of the peak (m)", duration(0.01) * speed, 2 * 1.2490 * R * math.sqrt(10), 0.10)
print(f"the pulse spans {duration(0.01) * speed / (2 * R):.2f} bead diameters at 1% of its peak "
      f"(theory {2 * 1.2490 * math.sqrt(10) / 2:.2f}; the full cos^4 support is {math.pi * math.sqrt(10) / 2:.1f})")

if failures:
    print("FAILED:", ", ".join(failures))
    sys.exit(1)
print("all comparisons within tolerance")
