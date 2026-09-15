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
derived from mesh sidesets and moving walls are a later layer.

The walls are given by the `wall_points` and `wall_normals` parameters of [KokkosParticleCloud.md],
one entry each per wall; normals are normalized on input.
