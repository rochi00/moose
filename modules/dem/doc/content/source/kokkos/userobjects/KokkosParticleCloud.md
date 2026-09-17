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
decision D5). Friction and rolling resistance are described in [Friction.md]; particle-wall
contacts take `wall_friction` and `wall_rolling_friction` when given, and the particle values
otherwise. The Coulomb limits of friction and rolling resistance are proportional to the magnitude of the net normal force (`friction_limit`, see [Friction.md]). With `limit_damping = true` the normal force is clamped at zero once the dashpot's
pull exceeds the spring's push (LAMMPS's option of the same name), so the contact never
attracts. With the Hertz model, `normal_damping` and `tangential_damping`, when given, replace
the Tsuji or Kuwabara-Kono damping by LAMMPS's `gran/hertz/history` dashpots
$\gamma\, m_\text{eff} \sqrt{r_\text{eff}\delta}\, v$. `rescale_histories` (default on) keeps
the length of a tangential or rolling spring rotated into the current tangent plane of a
turning contact, as LAMMPS's `pair granular` does, instead of only projecting it; and
`deformed_torque_arm` takes the torque of the tangential force at the contact plane,
$r_i - \delta/2$ from the center, as that style does, instead of at the radius. With
`mass_scaled_damping` the linear model's dashpot coefficients multiply the effective mass of the
contact, as the LAMMPS `gran/*` styles' $\gamma$ (a rate) does.

Particles listed in `initial_frozen` are never moved and act as immovable obstacles to the
others (LAMMPS `fix freeze`). With `disk_mass` the particles have the mass of disks, density
times $\pi r^2$, as LAMMPS gives them in two dimensions. Inserted particles take radii drawn
uniformly from `insertion_radius_range` when given (LAMMPS `fix pour diam range`). The
analytic walls can rotate rigidly about an axis (`wall_rotation_center`,
`wall_angular_velocity`): their points and normals turn with it and the surface velocity at
the contact point, which enters the contact, includes it.

The forces of the last substep stay with the particles between steps, migrating and restarting
with them, and drive the next step's first half-kick as they do between substeps; the forces
reported at the end of a step (by [KokkosParticleState.md] and the wall-force postprocessors)
are recomputed from the settled positions and the full-step velocities, so the trajectory does
not depend on how the substeps are grouped into steps.

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
its unwrapped position, just outside its rank's box, so its neighbors and their ghost images see
it where it is, until it is settled: wrapped by the period, placed by the point locator, and
migrated to the rank owning the far side, which need not be adjacent. Settling happens at the end
of every step and, within the substep loop, as soon as any particle is more than half the skin
outside its rank's box (across a periodic face or a partition), since beyond that its far-side
neighbors may not be ghosted here (this requires a replicated mesh, since the locator of a distributed mesh only sees the
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

## Insertion and outflow

With an `insertion_box` and `insertion_rate`, every step (between `insertion_start_time` and
`insertion_end_time`) the cloud inserts the step's share of the rate, the fraction being carried
over, at random positions in the box (uniform, in the mesh dimension's coordinates, from the
`insertion_seed`) that overlap no existing particle on any rank nor an earlier candidate of the
step, with the cloud's `radius` and `density`, `insertion_velocity`, and no spin, trying at most
ten candidates per particle. Every rank draws the same candidates and takes part in the
overlap decision, so the insertion is the same on any number of ranks; the rank whose element
holds an accepted position keeps it. New particles get the next global IDs. Particles whose
center lies in the `outflow_box` at the end of a step are removed, counted with those that left
the mesh in `num_exited`. The generator's position is checkpointed so a restart continues the
same sequence. The `pour` test inserts at 400 per second above a floor and drains through an
outflow box.

## Checkpoints and restart

The cloud is restartable: a checkpoint holds its local particles (as migration records, with the
libMesh ID of their element), its contact histories, the positions and velocities of moving
sideset walls, and its cumulative counters, filled from the device when the checkpoint is written.
On a restart or recovery the particles are placed back in their elements, the histories
reinserted, and the neighbor list, ghosts, and forces rebuilt, so the run continues as if
uninterrupted; a restep or backup restores the same way in place. The mesh, its partitioning, and
the number of ranks must be those of the checkpoint (a particle whose element is not on its rank
is an error). Particles unresolved at the checkpoint are located by the point locator.

## Checks

At setup the cloud errors on a mesh it cannot handle: a face that is not planar (a vertex more
than 1e-8 element sizes out of the plane of the face's first three vertices; the face walk is
exact only for planar faces). Mesh adaptivity is not supported. A particle with more
sideset wall contacts at once than can be evaluated (16) is an error rather than a loss of
contacts, and a particle left unresolved at the end of a step (not located in any element) is an
error unless `allow_unresolved = true`.

The substep is checked against the contact time scales at setup and every step, as LIGGGHTS's
`fix check/timestep/gran` does (`timestep_check`, `warn` by default, `error`, or `none`): for the
linear model the contact duration $\pi \sqrt{m_\text{eff} / k}$ of the stiffest spring, for the
Hertz model the Rayleigh time $\pi r_\text{min} \sqrt{\rho / G} / (0.1631 \nu + 0.8766)$ and the
Hertz contact time $2.87 (m_\text{eff}^2 / (r_\text{eff} E^{*2} v_\text{max}))^{1/5}$ at twice
the largest particle speed, with the smallest particle's mass and radius. The substep must be
below `timestep_fraction` (0.2 by default, plan Section 5.1) of each; the warning is given once.
Rolling springs are not included.

## Example Input Syntax

!listing modules/dem/test/tests/kokkos/integration/ballistic.i block=UserObjects

!syntax parameters /UserObjects/KokkosParticleCloud

!syntax inputs /UserObjects/KokkosParticleCloud

!syntax children /UserObjects/KokkosParticleCloud
