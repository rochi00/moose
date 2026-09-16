#!/usr/bin/env python3
"""Run oblique_impact.i over impact angles and plot the rebound spin and tangential restitution
against the gross-sliding solution, for the Hertz-Mindlin no-slip spring and the partial-slip
(Mindlin-Deresiewicz) spring. Usage: oblique_impact_sweep.py [dem-opt]"""
import csv, math, os, subprocess, sys

here = os.path.dirname(os.path.abspath(__file__))
exe = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "../../../../dem-opt")
v, r, mu, e_n = 3.85, 0.0025, 0.092, 0.98
angles = list(range(5, 90, 5))
theta_c = math.degrees(math.atan(3.5 * mu * (1 + e_n)))
print(f"critical angle for sliding throughout: {theta_c:.1f} degrees")

results = {}
for partial in (False, True):
    rows = []
    for deg in angles:
        th = math.radians(deg)
        vx, vy = v * math.sin(th), -v * math.cos(th)
        base = f"sweep_{'md' if partial else 'ns'}_{deg}"
        cmd = [exe, "-i", os.path.join(here, "oblique_impact.i"),
               f"UserObjects/cloud/initial_velocity='{vx} {vy} 0'",
               f"UserObjects/cloud/partial_slip={'true' if partial else 'false'}",
               f"Outputs/file_base={base}", "--error"]
        subprocess.run(cmd, cwd=here, check=True, capture_output=True)
        state = list(csv.reader(open(os.path.join(here, f"{base}_state_0003.csv"))))[1]
        vx1, vy1, wz1 = float(state[1]), float(state[2]), float(state[6])
        # Contact-point tangential velocities before and after: v_t + r omega (the spin about z
        # of a sphere on the floor below adds r omega_z to the x velocity of the contact point)
        e_t = -(vx1 + r * wz1) / (vx + 0)
        # Sliding-throughout prediction
        vn = -vy
        vx_s = vx - mu * (1 + e_n) * vn
        w_s = 2.5 * mu * (1 + e_n) * vn / r
        rows.append((deg, vx1, vy1, wz1, e_t, vx_s, w_s, vy1 / vn))
        for f in os.listdir(here):
            if f.startswith(base):
                os.remove(os.path.join(here, f))
    results[partial] = rows
    print(f"\n{'partial slip' if partial else 'no slip'}:")
    print(" angle   v_t'    omega'   e_t    | sliding v_t'  omega'    e_n")
    for row in rows:
        print(" %4d  %7.3f %8.1f  %6.3f |  %7.3f  %8.1f   %.4f" % (row[0], row[1], row[3], row[4], row[5], -row[6], row[7]))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for partial, style, label in ((False, "o-", "Hertz-Mindlin, no-slip spring"), (True, "s--", "partial slip (Mindlin-Deresiewicz)")):
        rows = results[partial]
        a1.plot([q[0] for q in rows], [q[3] for q in rows], style, label=label)
        a2.plot([q[0] for q in rows], [q[4] for q in rows], style, label=label)
    rows = results[False]
    # A sphere moving +x on the floor below spins about -z, and its contact point moves at
    # v_t + r w_z; sliding throughout gives v_t' + r w' = v_t - (7/2) mu (1 + e_n) v_n
    slide = [q for q in rows if q[0] > theta_c]
    a1.plot([q[0] for q in slide], [-q[6] for q in slide], "k:", label="sliding throughout (analytic)")
    a2.plot([q[0] for q in slide], [-1 + 3.5 * mu * (1 + e_n) / math.tan(math.radians(q[0])) for q in slide], "k:", label="sliding throughout (analytic)")
    for a in (a1, a2):
        a.axvline(theta_c, color="gray", lw=0.8)
        a.set_xlabel("impact angle from the normal (deg)")
        a.grid(alpha=0.3)
    a1.set_ylabel("rebound angular velocity (rad/s)")
    a2.set_ylabel("tangential restitution  -(v_t' + r w') / v_t")
    a2.set_ylim(-1.1, 1.1)
    a1.legend(fontsize=8)
    fig.suptitle("5 mm alumina sphere on glass at 3.85 m/s, mu = 0.092, e_n = 0.98 (Kharaz et al. 2001 setting)")
    fig.tight_layout()
    fig.savefig(os.path.join(here, "../../../../doc/content/media/dem_oblique_impact_sweep.png"), dpi=130)
    print("\nplot: doc/content/media/dem_oblique_impact_sweep.png")
except ImportError:
    pass
