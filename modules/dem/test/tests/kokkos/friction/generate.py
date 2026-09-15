#!/usr/bin/env python3
"""Generate the gold for sliding_sphere.i: a sphere launched sliding along a wall with Coulomb
friction decelerates at mu g while the friction torque spins it up, until the contact point stops
slipping at t* = 2 v0 / (7 mu g); it then rolls at 5/7 v0 (for a solid sphere, from the invariant
m v - I omega / r of the friction force acting at the contact point)."""
import math, os

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "gold"), exist_ok=True)

r, rho, g, mu, k = 0.01, 2500.0, 9.81, 0.3, 1e5
v0, x0 = 1.0, 0.1
dt, num_steps = 0.03, 10
T = dt * num_steps
m = rho * 4.0 / 3.0 * math.pi * r**3

t_roll = 2 * v0 / (7 * mu * g)
assert t_roll < T
v = 5.0 / 7.0 * v0
x = x0 + v0 * t_roll - 0.5 * mu * g * t_roll**2 + v * (T - t_roll)
# Rolling in +x on the floor below is a negative spin about z
omega = -v / r
# Rest on the linear spring: the center sits the static overlap m g / k below the radius
y = r - m * g / k
print(f"m = {m:.6e}, slipping ends at t = {t_roll:.6f} s; rolling at v = {v:.9f}, omega = {omega:.9f}")

with open(os.path.join(here, "gold", f"sliding_sphere_out_state_{num_steps:04d}.csv"), "w") as f:
    f.write("gid,vx,vy,vz,wx,wy,wz,x,y,z\n")
    f.write(f"0,{v:.17g},0,0,0,0,{omega:.17g},{x:.17g},{y:.17g},0\n")
