# SidesetWalls

Walls derived from mesh sidesets for a [KokkosParticleCloud.md] (plan layer L5, decision D9),
given by its `wall_boundaries` parameter. The faces of the listed boundaries on the local and
one-layer ghost elements are exported to device, triangulated (a quad is split along its first
diagonal; in 2D a face is an edge segment) with their normals pointing into the mesh, which is the
inside of the walls. Each local or ghost element lists the faces whose bounding box, inflated by
the largest radius plus the neighbor-list skin, meets its own, so a particle only tests the faces
near the element it is tracked in ([KokkosParticleCloud.md#element-tracking]). The lists are built
once at setup, the mesh being static. A particle whose element is unresolved sees no sideset
walls until it is placed at the end of the step.

## Narrow phase

For each candidate face the closest point to the sphere center is found (Ericson, *Real-Time
Collision Detection*, 5.1.5), giving a contact when it is closer than the radius. A contact with a
point in the face's interior uses the face normal and the signed distance to the face plane, so a
center that has passed behind the face is still pushed back through it; an edge or vertex contact
uses the direction from the point to the center. The contacts of a sphere are then reduced so that
a point shared by several faces gives one force:

- contacts with the same closest point are merged (the lowest-numbered face is kept), so a convex
  edge or vertex of the walls is a single contact along the line to the center;
- an edge or vertex contact is dropped when its point lies on a face that has an interior
  contact, so a sphere over a shared edge of a flat surface is not touched twice.

A triangulated plane thus acts exactly as a plane (the `plane3d` test slides a bouncing sphere
across triangle edges and a vertex shared by eight triangles and matches the analytic plane to
roundoff, plan case V12), and the reentrant corner of an L-shaped mesh gives the sphere-wall
restitution of a single contact (`corner`).

The contact acts as a rigid partner at rest of infinite radius and mass, with the same normal,
tangential, and rolling models as sphere-sphere contact ([LinearSpringDashpot.md], [Hertz.md],
[Friction.md]). The contact history is keyed by the face's boundary, not the face, so it carries
across the faces of one sideset as the contact point slides over them; faces a particle can touch
at the same time with different normals (the floor and a side of a box) should therefore be in
different sidesets. A boundary's history is advanced once per substep, from its first overlapping
face, and reset only when none of its listed faces overlaps.
