# SidesetWalls

Walls derived from mesh sidesets for a [KokkosParticleCloud.md] (plan layer L5, decision D9),
given by its `wall_boundaries` parameter. The faces of the listed boundaries are triangulated (a
quad is split along its first diagonal; in 2D a face is an edge segment) with their normals
pointing into the mesh, which is the inside of the walls, and exported to device. Every rank
holds the faces of its own elements and, as LIGGGHTS distributes its meshes, receives from the
other ranks the faces within a particle's reach (the largest radius plus the neighbor-list skin
and `wall_reach`) of its bounding box, so a particle can touch any face it can reach whatever the
partitioning. Each local or ghost element lists the faces whose bounding box, inflated by the
reach, meets its own, so a particle only tests the faces near the element it is tracked in
([KokkosParticleCloud.md#element-tracking]). The lists are built once at setup, the mesh being
static. A particle whose element is unresolved sees no sideset walls until it is placed at the
end of the step.

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
[Friction.md]). The contact history is keyed by the face. Faces sharing an edge (a node in 2D)
whose normals are within `wall_curvature` degrees of each other belong to one surface, and when
the contact point slides from one to the other the history is handed over (LIGGGHTS's
`curvature`): a touching face with an unstretched history takes the stretched history of a
same-surface neighbor the contact has just left. So a triangulated plane keeps its tangential
spring like a plane, while faces meeting at a crease keep separate histories, and the sides of a
box may share a sideset (the `hopper_one_sideset` test). Histories are kept for every face within
the skin, listed without the contact reduction so that a face whose edge contact is hidden by a
neighbor's face contact, the next face to be touched, already has one.

## Moving walls

With `wall_displacements`, nodal displacement variables (one per mesh dimension, for example
driven by a `FunctionAux` or a solid mechanics solve), the walls move: at every MOOSE step the
faces are carried from where they are to the displaced positions of their nodes at the constant
velocity that gets them there over the step, advancing every substep along with the particles,
and the velocity of a contact point, interpolated from the face's vertices, enters the relative
velocity of the contact (its dashpots and friction) like a partner's. The face normals are
recomputed for the new positions. The candidate lists of the elements are not rebuilt, so the
motion must stay within `wall_reach` of the undisplaced faces, and the particles are still tracked
on the undisplaced mesh. The `moving_floor` test (plan case V14) rests a sphere on a floor driven
sideways and upward and checks that it rides along to 1e-8.
