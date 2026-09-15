# KokkosParticleCloud

!syntax description /UserObjects/KokkosParticleCloud

## Overview

`KokkosParticleCloud` holds a cloud of spherical particles on device and, every time it executes,
integrates them over `substeps` substeps of the MOOSE time step with Velocity Verlet (translation)
and an exact constant-angular-velocity quaternion rotation (orientation), then re-resolves the
local element containing each particle. Particle state is stored in Kokkos views using the libMesh
Kokkos storage policy (`libmesh/kokkos_storage_policy.h`). The substep loop runs entirely on
device (plan decision D4).

Forces are the body force `m * gravity` and the normal contact force of the selected
`contact_model`, [LinearSpringDashpot.md] (the default, enabled by a positive `normal_stiffness`)
or [Hertz.md], between overlapping spheres found through the neighbor list, against the fixed
planar walls given by `wall_points` and `wall_normals` ([AnalyticWalls.md]), and against the
sidesets listed in `wall_boundaries` ([SidesetWalls.md]). The model is a compile-time policy: the
force kernel is instantiated once per model and the selection is made once at setup (plan
decision D5). Friction and rolling resistance are described in [Friction.md].

## Neighbor list

`DEM::NeighborList` lists every pair closer than `r_i + r_j + skin` from both sides, so force
kernels can sum a particle's contacts without atomics. The candidate pairs come from the
`broad_phase` (plan decision D7):

- `uniform_grid` bins the particles on a grid of cell size `2 r_max + skin` and searches the 27
  surrounding cells. It is the fastest for a narrow size distribution, but the cell size is set
  by the largest particle, so it degrades as the size ratio grows.
- `arborx_bvh` builds an [ArborX](https://github.com/arborx/ArborX) bounding volume hierarchy
  over the spheres' bounding boxes and queries it with each box inflated by the skin; the cost
  does not depend on the size distribution. It needs the ArborX submodule
  (`modules/dem/contrib/arborx`, header-only; it requires C++20, so the module's Kokkos sources
  are compiled as C++20 when it is present) and is reported by the `arborx` capability.

Each particle's neighbors are then sorted by index, so the list, and with it the summation order
of the forces, is the same whichever broad phase built it. The list is rebuilt whenever, on any
rank, some particle has moved more than `skin / 2` since the last build or the number of local
particles changed. With `verify_neighbor_list = true`, every build is checked against a
brute-force pair search on host.

## Ghost particles

Contact across partitions uses ghost particles: at every rebuild, each rank sends copies of its
particles lying within the neighbor-list cutoff `2 r_max + skin` of a neighboring rank's inflated
bounding box, and appends the copies it receives after its own particles. The neighbor list covers
local and ghost particles, and every rank computes the forces on its own particles from the ghost
copies, so no force communication is needed (a pair spanning two ranks is evaluated once on each).
Between rebuilds only the ghost positions and velocities are forwarded, every substep, with
nonblocking MPI: the forces on the local particles without ghost neighbors are computed while the
forward is in flight, and those with ghost neighbors once it has landed. The device buffers are
passed to MPI directly when the host can access device memory or `gpu_aware_mpi = true`, and
through host mirrors otherwise. A pair with
a ghost is counted on the rank owning the smaller global ID, so `num_neighbor_pairs` is the same on
any number of ranks; it is the count at the last build, which in parallel can be later in the step
than in serial because migration forces a rebuild.

## Pair states

Contacts carry a history, the tangential and rolling spring displacements of [Friction.md], in a
device hash
map (`Kokkos::UnorderedMap`, plan layer L3, decision D6) keyed by the global IDs of the pair in
increasing order, or by the particle's ID and the wall's index for a wall contact. The history is
stored in the lower-ID particle's orientation, so every rank holding either particle forms the
same key and, advancing it from the same ghost positions, velocities, and spins, keeps an
identical copy: no history or force communication is needed for a pair across a partition.

Every substep, each pair is advanced from exactly one side per rank (a pair of local particles
from the lower index, a pair with another rank's ghost from the local particle, a pair with a
periodic image of a local particle from the lower ID), and a listed pair that no longer overlaps
has its history reset; then every particle sums its forces from both sides of its contacts. The
map is rebuilt with the neighbor list, keeping the histories of the pairs still listed and the
wall contacts within the skin and dropping the rest, and a migrating particle carries the
stretched histories of its contacts to its new rank. With `verify_pair_states = true`, every step
checks that the histories with a stretched spring are exactly those of overlapping pairs and that
every listed pair has a history. `num_pair_states` of [KokkosParticleCloudValue.md] counts the
histories held, summed over ranks.

The `potential_energy` reported excludes the tangential and rolling springs.

## Periodic directions

With `periodic`, the domain wraps over the extent of the mesh bounding box in the listed
directions. A particle that walks out of the mesh is then unresolved rather than exited: it keeps
its unwrapped position, just outside its rank's box, for the rest of the step, so its neighbors and
their ghost images see it where it is, and at the end of the step it is wrapped by the period,
placed by the point locator, and migrated to the rank owning the far side, which need not be
adjacent (this requires a replicated mesh, since the locator of a distributed mesh only sees the
local and ghost elements). Contact across a periodic face uses the ghost machinery: every rank
sends, to every rank including itself, the particles whose image under each combination of period
shifts lies within the cutoff of that rank's box, with the shifted positions, so the neighbor list
needs no minimum image convention. A rank's box must be wider than the cutoff so a particle does
not interact with its own image.

## Element tracking

Element tracking is incremental: each particle remembers its element, and after moving it is tested
against the outward face planes cached by the Kokkos mesh (`Moose::Kokkos::Mesh::getSideNormal`,
`getSideCentroid`). A particle in front of a face plane leaves through the face it is furthest in
front of and the test repeats on the face neighbor, up to `max_hops` times per substep. This is
exact only for elements with planar faces.

The walk ends in one of four states:

- inside a local element;
- inside a one-layer ghost element, which has no cached face geometry. The particle keeps
  integrating with that stale element until the end of the MOOSE step, when it is sent to the rank
  owning the element, which resumes the walk from there. Rounds of walking and exchanging repeat
  until no particle is left in a ghost element, so a particle may cross several partitions in one
  step. Departing particles are found, compacted out of the cloud, and packed on device, so only
  the departing and arriving particles cross to the host for the TIMPI exchange; GPU-aware MPI is
  deferred;
- exited, because it walked out of the mesh. An exited particle is inert for the rest of the step
  and removed at the end of it;
- unresolved, because it exhausted `max_hops`. An unresolved particle keeps integrating and is
  located again with libMesh's `PointLocator` at the end of the step.

With `verify = true`, each step the device assignment is checked on host against libMesh's
`PointLocator` and any mismatch is an error.

Global scalars are exposed through [KokkosParticleCloudValue.md] and per-particle state through
[KokkosParticleState.md].

## Example Input Syntax

!listing modules/dem/test/tests/kokkos/integration/ballistic.i block=UserObjects

!syntax parameters /UserObjects/KokkosParticleCloud

!syntax inputs /UserObjects/KokkosParticleCloud

!syntax children /UserObjects/KokkosParticleCloud
