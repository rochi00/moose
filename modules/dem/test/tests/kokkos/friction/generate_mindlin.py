#!/usr/bin/env python3
"""Generate the gold for mindlin_loading.i (plan case V7): a sphere resting on a wall under gravity
with the Hertz-Mindlin model and partial slip is given a small tangential velocity. With the
normal force P = m g constant, the tangential force follows the Mindlin-Deresiewicz loading curve
  Q(delta) = mu P [1 - (1 - delta / delta_max)^(3/2)],  delta_max = 3 mu P / (16 G* a),
a = sqrt(r delta_n) being the contact radius, and, elastic on unloading, the sphere oscillates on
it. The motion is that of the contact point: its displacement delta obeys
(2/7) m delta'' = -Q(delta) (the force at the contact point drives the translation and, through
the torque r Q with I = (2/5) m r^2, the spin, so the contact point sees 1 + 5/2 times the
acceleration), the velocity is v0 - J / m and the spin -(5/2) J / (m r) for the impulse J. The
reference is integrated with fourth-order Runge-Kutta at a step where its error is far below
the tolerance. The launch speed is chosen so the excursion reaches 0.7 delta_max, well into
partial slip but short of gross sliding."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

r, rho, g, E, nu, mu = 0.01, 2500.0, 9.81, 1e7, 0.3, 0.3
x0 = 0.05
dt, num_steps = 0.002, 10
m = rho * 4.0 / 3.0 * math.pi * r**3
E_star = E / (2 * (1 - nu**2))
G_star = E / (2 * (1 + nu)) / (2 * (2 - nu))
P = m * g
delta_n = (P / (4.0 / 3.0 * E_star * math.sqrt(r)))**(2.0 / 3.0)
a = math.sqrt(r * delta_n)
k_t = 8 * G_star * a
delta_max = 3 * mu * P / (2 * k_t)
y = r - delta_n

def Q(d):
    s = min(abs(d), delta_max) / delta_max
    return math.copysign(mu * P * (1 - (1 - s)**1.5), d)

# Launch speed for a 0.7 delta_max excursion: the contact-point kinetic energy (2/7) m u^2 / 2
# equals the work of Q up to there
target = 0.7 * delta_max
work = mu * P * delta_max * (0.7 - 0.4 * (1 - (1 - 0.7)**2.5))
u0 = math.sqrt(2 * work / (2.0 / 7.0 * m))
v0 = u0                                  # no initial spin, so the contact point moves with the center
print(f"m = {m:.6e}, delta_n = {delta_n:.6e}, a = {a:.6e}, k_t = {k_t:.6f}, "
      f"delta_max = {delta_max:.6e}, v0 = {v0!r}, period ~ {2 * math.pi * math.sqrt(2 / 7 * m / k_t):.6f} s")

# RK4 on (delta, u, J, x): delta' = u, u' = -Q(delta) / (2 m / 7), J' = Q(delta), x' = v0 - J / m
m_c = 2.0 / 7.0 * m
def rhs4(s):
    d, u, J, x = s
    return (u, -Q(d) / m_c, Q(d), v0 - J / m)
h = 1e-7
t = 0.0
state = (0.0, u0, 0.0, x0)
golds = {}
for step in range(1, num_steps + 1):
    while t < step * dt - 0.5 * h:
        k1 = rhs4(state)
        k2 = rhs4(tuple(s + 0.5 * h * k for s, k in zip(state, k1)))
        k3 = rhs4(tuple(s + 0.5 * h * k for s, k in zip(state, k2)))
        k4 = rhs4(tuple(s + h * k for s, k in zip(state, k3)))
        state = tuple(s + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4) for s, a1, a2, a3, a4 in zip(state, k1, k2, k3, k4))
        t += h
    d, u, J, x = state
    v = v0 - J / m
    omega = -2.5 * J / (m * r)
    print(f"t = {step * dt:.3f}: delta / delta_max = {d / delta_max:+.4f}, Q / mu P = {Q(d) / (mu * P):+.4f}")
    golds[step] = (v, omega, x, -Q(d))
for step in (5, 10):
    v, omega, x, fx = golds[step]
    with open(os.path.join(here, "gold", f"mindlin_loading_out_state_{step:04d}.csv"), "w") as f:
        # The torque of the tangential force at the contact point below the center is r fx about z
        f.write("fx,fy,fz,gid,taux,tauy,tauz,vx,vy,vz,wx,wy,wz,x,y,z\n")
        f.write(f"{fx:.17g},0,0,0,0,0,{r * fx:.17g},{v:.17g},0,0,0,0,{omega:.17g},{x:.17g},{y:.17g},0\n")
print(f"initial position = {x0!r} {y!r}")
