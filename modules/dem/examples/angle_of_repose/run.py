#!/usr/bin/env python3
"""Angle of repose of coarse spheres (Zhou, Xu, Yu & Zulli 2002, Powder Technol. 125, 45-54):
a bed of 10 mm spheres settled on a 30 d plate in a container 40 d wide and w thick is
discharged over both edges of the plate, leaving a heap whose surface slope is the angle of
repose. Zhou et al. fit their DEM results (standard deviation 2 degrees) with

    theta  = 102   mu_s,pp^0.27 mu_s,pw^0.22 mu_r,pp^0.06 mu_r,pw^0.12 d^-0.2   (w = 4 d, Eq. 3)
    theta0 = 68.61 mu_s,pp^0.27 mu_s,pw^0.22 mu_r,pp^0.06 mu_r,pw^0.12 d^-0.2   (w -> inf, Eq. 5)

with d in mm and their rolling friction coefficients in mm (torque mu_r F_n); ours are
dimensionless with torque mu_r r_eff F_n, r_eff = r/2 between particles and r against walls, so
mu_r,pp = 0.05 mm and mu_r,pw = 0.1 mm both map to mu_r = 0.02 for d = 10 mm. Their experiments
with glass beads give theta0 = 29.6 d^-0.206 (18.4 degrees at 10 mm) and are matched by
mu_s,pp = 0.4, mu_r,pp = 0.05 mm with the wall values 1.5 and 2 times those.

Each case settles a random rain of spheres on the closed plate (stage 1), then transplants the
bed onto the plate with open outlets (stage 2) and reads the slope of both flanks of the heap.
Usage: run.py [-n ranks] [dem-opt] [case names...]"""
import csv, glob, math, os, subprocess, sys

here = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]
ranks = 8
if args[:1] == ["-n"]:
    ranks = int(args[1])
    args = args[2:]
exe = args[0] if args and os.path.sep in args[0] else os.path.join(here, "../../dem-opt")
names = [a for a in args if os.path.sep not in a]

d = 0.010
r = d / 2
rho = 2500
g = 9.81
m = rho * 4 / 3 * math.pi * r**3
# Linear spring: a 15 d column overlaps the bottom spheres by 15 m g / k = 0.4% of r; damping
# gives a restitution of about 0.45 (Zhou et al.: the angle is insensitive to stiffness and
# damping); the substep is a 27th of the contact time pi sqrt(m_eff / k), m_eff = m / 2
k = 1e4
c = 0.25 * 2 * math.sqrt(k * m / 2)
substeps = 333
plate = 15 * d                    # half-width of the plate (30 d)
box = 20 * d                      # half-width of the container (40 d)
height = 20 * d

# name: (thickness / d, periodic in y, mu_pp, mu_pw, mu_r_pp, mu_r_pw, particles, Zhou reference)
def zhou(mu_s_pp, mu_s_pw, mu_r_pp, mu_r_pw, thick):
    """Eq. 3 (w = 4 d) or Eq. 5 (w -> inf) with the rolling coefficients converted to mm"""
    A = 102 if thick == 4 else 68.61
    return A * mu_s_pp**0.27 * mu_s_pw**0.22 * (mu_r_pp * r / 2 * 1e3)**0.06 * (mu_r_pw * r * 1e3)**0.12 * (d * 1e3)**-0.2

# name: (thickness / d, periodic in y, mu_pp, mu_pw, mu_r_pp, mu_r_pw, particles). Sliding
# friction as Zhou et al.'s base condition; the rolling coefficient is swept because their
# rolling model (a constant torque against each particle's own spin) is not the EPSD spring
# used here, so their mu_r does not carry over; the thickness series then tests the
# experimental wall effect theta = theta0 (1 + exp(-0.18 w/d)) at the mu_r that matches the
# 4 d experiment (27.4 degrees for 10 mm glass beads)
cases = {}
for mu_r in (0.02, 0.05, 0.1, 0.2):
    cases[f"w4_mur{mu_r:g}"] = (4, False, 0.4, 0.6, mu_r, mu_r, 2000)
for mu_r in (0.05, 0.1):
    cases[f"w8_mur{mu_r:g}"] = (8, False, 0.4, 0.6, mu_r, mu_r, 4000)
    cases[f"w8p_mur{mu_r:g}"] = (8, True, 0.4, 0.6, mu_r, mu_r, 4000)
names = names or list(cases)


def run(base):
    cmd = ["mpiexec", "-n", str(ranks), exe, "-i", base + ".i"] if ranks > 1 else [exe, "-i", base + ".i"]
    out = subprocess.run(cmd, cwd=here, capture_output=True, text=True)
    if out.returncode:
        print(out.stdout[-3000:], out.stderr[-2000:])
        sys.exit(1)


def cloud(positions, walls, mu, mu_w, mu_r, mu_r_w, periodic, extra=""):
    return f"""[UserObjects]
  [cloud]
    type = KokkosParticleCloud
    initial_positions = '{chr(10).join(f"{x:.6f} {y:.6f} {z:.6f}" for x, y, z in positions)}'
    radius = {r}
    density = {rho}
    gravity = '0 0 -{g}'
    normal_stiffness = {k}
    normal_damping = {c}
    tangential_stiffness = {k}
    tangential_damping = {c}
    friction = {mu}
    wall_friction = {mu_w}
    rolling_stiffness = {k}
    rolling_damping = {c}
    rolling_friction = {mu_r}
    wall_rolling_friction = {mu_r_w}
    wall_boundaries = '{walls}'
{"    periodic = y" + chr(10) if periodic else ""}    skin = {r / 2}
    substeps = {substeps}
    timestep_check = error
    execute_on = TIMESTEP_END
{extra}  []
[]
[VectorPostprocessors]
  [state]
    type = KokkosParticleState
    cloud = cloud
    execute_on = TIMESTEP_END
  []
[]
[Postprocessors]
  [ke]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = kinetic_energy
  []
  [num]
    type = KokkosParticleCloudValue
    cloud = cloud
    value = num_particles
  []
[]
"""


def mesh(half_width, thick, periodic):
    return f"""[Mesh]
  [box]
    type = GeneratedMeshGenerator
    dim = 3
    nx = {int(round(2 * half_width / (2 * d)))}
    ny = {max(2, int(round(thick / (2 * d))))}
    nz = {int(round(height / (2 * d)))}
    xmin = {-half_width:.4f}
    xmax = {half_width:.4f}
    ymax = {thick:.4f}
    zmax = {height:.4f}
  []
"""


outputs = """[Outputs]
  [csv]
    type = CSV
    execute_postprocessors_on = TIMESTEP_END
    execute_vector_postprocessors_on = FINAL
  []
  [console]
    type = Console
    execute_postprocessors_on = NONE
  []
[]
"""

results = []
for name in names:
    wd, periodic, mu, mu_w, mu_r, mu_r_w, N = cases[name]
    thick = wd * d
    # GeneratedMesh names: z = 0 is 'back', the y faces are 'bottom' and 'top', the x faces 'left' and 'right'
    walls = "back left right" if periodic else "back left right bottom top"
    # Stage 1: rain N spheres onto the closed plate and let them settle
    rate = 4000
    t1 = N / rate + 0.8
    with open(os.path.join(here, name + "_fill.i"), "w") as f:
        f.write(mesh(plate, thick, periodic) + "[]\n[Problem]\n  solve = false\n[]\n" +
                cloud([(0, thick / 2, r * 1.05)], walls, mu, mu_w, mu_r, mu_r_w, periodic, f"""    insertion_box = '{-plate + 1.1 * r:.5f} {1.1 * r if not periodic else 0:.5f} {height - 4 * d:.5f}
                     {plate - 1.1 * r:.5f} {thick - 1.1 * r if not periodic else thick:.5f} {height - 1.1 * r:.5f}'
    insertion_rate = {rate}
    insertion_velocity = '0 0 -1'
    insertion_end_time = {N / rate:.4f}
    insertion_seed = 11
""") + f"""[Executioner]
  type = Transient
  dt = 0.01
  num_steps = {int(round(t1 / 0.01))}
[]
""" + outputs)
    if not glob.glob(os.path.join(here, name + "_fill_csv_state_*.csv")):
        print(f"{name}: settling {N} spheres on the closed plate ...", flush=True)
        run(name + "_fill")
    state = sorted(glob.glob(os.path.join(here, name + "_fill_csv_state_*.csv")))[-1]
    bed = [(float(q["x"]), float(q["y"]), float(q["z"])) for q in csv.DictReader(open(state))]
    # Stage 2: the same bed on the plate inside the wider container; the bottom faces beyond
    # the plate are left open (not walls), so spheres discharging over the edges leave the mesh
    with open(os.path.join(here, name + ".i"), "w") as f:
        f.write(mesh(box, thick, periodic) + f"""  [plate]
    type = ParsedGenerateSideset
    input = box
    combinatorial_geometry = 'abs(x) < {plate:.5f}'
    included_boundaries = back
    new_sideset_name = plate
    replace = true
  []
[]
[Problem]
  solve = false
[]
""" + cloud(bed, walls.replace("back", "plate"), mu, mu_w, mu_r, mu_r_w, periodic) + f"""[Executioner]
  type = Transient
  dt = 0.01
  num_steps = 600
[]
""" + outputs)
    run(name)
    final = sorted(glob.glob(os.path.join(here, name + "_csv_state_*.csv")))[-1]
    heap = [(float(q["x"]), float(q["z"])) for q in csv.DictReader(open(final))]
    hist = list(csv.DictReader(open(os.path.join(here, name + "_csv.csv"))))
    ke = [float(q["ke"]) for q in hist]
    # Surface profile: the top sphere in each column of width d; slopes fitted on both flanks
    # between 15% and 85% of the plate half-width, as read off Zhou et al.'s profiles
    nb = int(2 * plate / d)
    top = [None] * nb
    for x, z in heap:
        b = min(max(int((x + plate) / d), 0), nb - 1)
        top[b] = z if top[b] is None else max(top[b], z)
    slopes = []
    for side in (-1, 1):
        pts = [((b + 0.5) * d - plate, top[b]) for b in range(nb)
               if top[b] is not None and 0.15 * plate < side * ((b + 0.5) * d - plate) < 0.85 * plate]
        xm, zm = sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
        s = sum((p[0] - xm) * (p[1] - zm) for p in pts) / sum((p[0] - xm) ** 2 for p in pts)
        slopes.append(math.degrees(math.atan(abs(s))))
    theta = sum(slopes) / 2
    # Zhou et al.'s glass-bead experiments: theta0 = 29.6 d^-0.206 without front and rear walls,
    # times 1 + exp(-0.18 w/d) with them (their Eqs. 4 and 13 fit)
    theta0_exp = 29.6 * (d * 1e3) ** -0.206
    exp = theta0_exp * (1 if periodic else 1 + math.exp(-0.18 * wd))
    ref = zhou(mu, mu_w, mu_r, mu_r_w, wd if not periodic else 0)
    results.append((name, wd, periodic, mu_r, theta, slopes, exp, ref, len(heap)))
    print(f"{name}: angle of repose {theta:.1f} deg (flanks {slopes[0]:.1f}, {slopes[1]:.1f}); glass-bead "
          f"experiment {exp:.1f} deg ({theta - exp:+.1f}); Zhou et al.'s DEM fit with their rolling model "
          f"{ref:.1f}; {len(heap)} of {len(bed)} spheres left on the plate, KE {ke[-1]:.1e} J (peak {max(ke):.1e})",
          flush=True)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(9, 4))
    w4 = sorted((q[3], q[4]) for q in results if q[1] == 4 and not q[2])
    axs[0].plot([q[0] for q in w4], [q[1] for q in w4], "o-", label="DEM, w = 4 d")
    axs[0].axhline(29.6 * 10 ** -0.206 * (1 + math.exp(-0.72)), color="k", ls="--",
                   label="10 mm glass beads, w = 4 d (Zhou et al. 2002)")
    axs[0].set_xlabel("rolling friction coefficient mu_r")
    axs[0].set_ylabel("angle of repose (deg)")
    for mu_r, mk in ((0.05, "o"), (0.1, "s")):
        pts = sorted((q[1] if not q[2] else 40, q[4]) for q in results if q[3] == mu_r)
        if pts:
            axs[1].plot([p[0] for p in pts], [p[1] for p in pts], mk, label=f"DEM, mu_r = {mu_r:g} (periodic drawn at 40)")
    ww = [2 + 0.5 * i for i in range(77)]
    axs[1].plot(ww, [29.6 * 10 ** -0.206 * (1 + math.exp(-0.18 * w)) for w in ww], "k--",
                label="experiment: theta0 (1 + exp(-0.18 w/d))")
    axs[1].set_xlabel("container thickness w / d")
    axs[1].set_ylabel("angle of repose (deg)")
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(here, "angle_of_repose.png"), dpi=130)
    print("plot: angle_of_repose.png")
except ImportError:
    pass
