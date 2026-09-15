#!/usr/bin/env python3
"""Generate the golds for incline.i and rolling_stop.i (plan case V8).

incline.i: a sphere placed at rest on the plane y = 0.1 - x tan(theta), theta = 10 degrees,
already pressed in by the static normal overlap so the friction holds from the start.
  - incline_rest: mu = 0.5 and mu_r = 0.25, both above tan(theta) = 0.176. The friction holds
    the sphere and the rolling resistance torque mu_r m g cos(theta) r can balance the torque
    r m g sin(theta) of the friction, so it stays at rest, settled into the tangential and
    rolling springs (each stretched by m g sin(theta) over its stiffness; a rolling stretch moves
    the center by as much).
  - incline_roll: mu_r = 0.1 < tan(theta) < mu. The rolling spring saturates and the sphere
    rolls without slipping under the plastic torque, at a = g (sin - mu_r cos) / (1 + 2/5).
  - incline_slide: mu = 0.15 < tan(theta) < mu_r = 0.25. The rolling resistance locks the
    rotation and the friction saturates: the sphere slides without rolling at
    a = g (sin - mu cos), the critical angle of a sliding body being atan(mu).
The golds of the two moving cases are exact for the steady acceleration; the start-up, in which
the sphere accelerates at g sin(theta) while the springs stretch from zero, leaves a velocity
offset of 3e-4 m/s that the test tolerance covers.

rolling_stop.i: a sphere rolling without slipping on a horizontal plane at v0, with mu_r = 0.1.
The rolling spring saturates at once (the contact point is at rest, so there is no tangential
transient) and the plastic torque decelerates it at mu_r g / (1 + 2/5) until it stops at
t = 1.4 v0 / (mu_r g) after v0^2 / (2 a), where the spring then holds it."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

r, rho, g, k_n, k_t, k_r = 0.01, 2500.0, 9.81, 1e5, 1e5, 1e5
m = rho * 4.0 / 3.0 * math.pi * r**3


def write(name, num_steps, v, omega, x):
    with open(os.path.join(here, "gold", f"{name}_out_state_{num_steps:04d}.csv"), "w") as f:
        f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
        f.write(f"0,{v[0]:.17g},{v[1]:.17g},0,0,0,{omega:.17g},{x[0]:.17g},{x[1]:.17g},0\n")


# --- incline
theta = math.radians(10)
contact0 = (0.05, 0.1 - 0.05 * math.tan(theta))
n = (math.sin(theta), math.cos(theta))         # plane normal, into the domain
t = (math.cos(theta), -math.sin(theta))        # downslope tangent
dt, num_steps = 0.05, 10
T = dt * num_steps
d_n = m * g * math.cos(theta) / k_n
x0 = [contact0[c] + (r - d_n) * n[c] for c in range(2)]
print(f"m = {m:.6e}, tan(theta) = {math.tan(theta):.6f}, normal = {n[0]!r} {n[1]!r}, "
      f"initial position = {x0[0]!r} {x0[1]!r}")
for name, mu, mu_r in (("incline_rest", 0.5, 0.25), ("incline_roll", 0.5, 0.1), ("incline_slide", 0.15, 0.25)):
    if mu_r < math.tan(theta):                 # rolls
        a = g * (math.sin(theta) - mu_r * math.cos(theta)) / 1.4
        rolling = True
    elif mu < math.tan(theta):                 # slides
        a = g * (math.sin(theta) - mu * math.cos(theta))
        rolling = False
    else:
        a = 0
        rolling = False
    s = 0.5 * a * T**2
    v = a * T
    # Quasi-static spring stretches: the tangential spring carries the friction force, the
    # rolling spring the friction torque, and each moves the center along the slope
    f_t = m * (g * math.sin(theta) - a)
    stretch = f_t / k_t + (f_t / k_r if not rolling and mu_r >= math.tan(theta) else 0)
    x = [x0[c] + (s + stretch) * t[c] for c in range(2)]
    omega = -v / r if rolling else 0
    print(f"{name}: a = {a:.9f}, v = {v:.9f}, omega = {omega:.9f}")
    write(name, num_steps, (v * t[0], v * t[1]), omega, x)

# --- rolling stop
v0, mu_r, x0 = 1.0, 0.1, 0.1
dt, num_steps = 0.1, 20
T = dt * num_steps
a = mu_r * g / 1.4
t_stop = v0 / a
assert t_stop < T
print(f"rolling_stop: a = {a:.9f}, stops at t = {t_stop:.6f} s after {v0**2 / (2 * a):.9f} m")
write("rolling_stop", num_steps, (0, 0), 0, (x0 + v0**2 / (2 * a), r - m * g / k_n))
