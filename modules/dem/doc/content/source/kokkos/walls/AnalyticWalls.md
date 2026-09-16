# AnalyticWalls

Fixed planar walls of a [KokkosParticleCloud.md] given analytically (plan layer L5), each by a
point on it and its normal pointing into the domain. The sphere-plane narrow phase is the signed
distance of the particle center to the plane: a sphere of radius $r$ at $x$ overlaps wall $w$ by

!equation
\delta = r - (x - p_w) \cdot n_w,

and while $\delta > 0$ it gets the [LinearSpringDashpot.md] normal force with $v_n = v \cdot n_w$,
since the wall is rigid and at rest. The wall therefore acts as a contact partner of infinite mass,
so the effective mass of a sphere-wall contact is the sphere's own and its coefficient of
restitution is that of the damped oscillator with $m_\text{eff} = m$, which the `wall_restitution`
test checks for an oblique impact (the tangential velocity is unchanged, the contact being
frictionless).

A wall extends without bound and is not tied to the mesh: the mesh must keep the particle centers
on the inner side of every wall, or a particle exits the mesh before it reaches the wall. Walls
tied to the mesh are [SidesetWalls.md].

The walls are given by the `wall_points` and `wall_normals` parameters of [KokkosParticleCloud.md],
one entry each per wall; normals are normalized on input.

## Moving and servo-controlled walls

A wall may translate at a prescribed `wall_velocities` entry, which its contacts see as the
partner's velocity and which carries its point along every substep (velocity-controlled shear,
compression). With `servo_walls`, `servo_forces`, `servo_gain`, and `servo_max_velocity`, a wall
is instead driven to a target compressive force, the particles' push against its inward normal
as measured by [KokkosWallForce.md] at the end of the previous step: every step its velocity
along the normal is the gain times the force error, capped by `servo_max_velocity` and by the
velocity that would change a single contact's force by the error within the step (from the
contact stiffness at the target force), which keeps the control stable for any gain and is what
limits the approach speed of a wall not yet in contact. The `compress` test squeezes a sphere
between the floor and a servo wall to 1 N, reached to 4e-11. Wall positions are checkpointed.
