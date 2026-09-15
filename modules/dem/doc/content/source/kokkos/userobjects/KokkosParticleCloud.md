# KokkosParticleCloud

!syntax description /UserObjects/KokkosParticleCloud

## Overview

`KokkosParticleCloud` holds a cloud of spherical particles on device and, every time it executes,
integrates them over `substeps` substeps of the MOOSE time step with Velocity Verlet (translation)
and an exact constant-angular-velocity quaternion rotation (orientation), then re-resolves the
local element containing each particle. Particle state is stored in Kokkos views using the libMesh
Kokkos storage policy (`libmesh/kokkos_storage_policy.h`). The substep loop runs entirely on
device (plan decision D4).

Forces are the body force `m * gravity` and, when `normal_stiffness` is positive, the
[LinearSpringDashpot.md] normal contact force between overlapping spheres found through the
neighbor list. Tangential contact, rolling resistance, and walls are added by later layers.

## Neighbor list

A uniform-grid broad phase (`DEM::NeighborList`) bins the local particles on a grid of cell size
`2 r_max + skin`, counting-sorts them by bin with a unique key so the order is deterministic, and
lists every pair closer than `r_i + r_j + skin` from both sides, so force kernels can sum a
particle's contacts without atomics. The list is rebuilt whenever some particle has moved more than
`skin / 2` since the last build or the number of local particles changed. Pairs spanning a
partition boundary are not listed until ghost particles are exchanged. With
`verify_neighbor_list = true`, every build is checked against a brute-force pair search on host.

## Element tracking

Element tracking is incremental: each particle remembers its element, and after moving it is tested
against the outward face planes cached by the Kokkos mesh (`Moose::Kokkos::Mesh::getSideNormal`,
`getSideCentroid`). A particle in front of a face plane leaves through the face it is furthest in
front of and the test repeats on the face neighbor, up to `max_hops` times per substep. This is
exact only for elements with planar faces.

The walk ends in one of three states:

- inside a local element;
- inside a one-layer ghost element, which has no cached face geometry. The particle keeps
  integrating with that stale element until the end of the MOOSE step, when it is sent to the rank
  owning the element, which resumes the walk from there. Rounds of walking and exchanging repeat
  until no particle is left in a ghost element, so a particle may cross several partitions in one
  step. Departing particles are found, compacted out of the cloud, and packed on device, so only
  the departing and arriving particles cross to the host for the TIMPI exchange; GPU-aware MPI is
  deferred;
- lost, because it walked out of the mesh or exhausted `max_hops`. Lost particles keep integrating.

With `verify = true`, each step the device assignment is checked on host against libMesh's
`PointLocator` and any mismatch is an error.

Global scalars are exposed through [KokkosParticleCloudValue.md] and per-particle state through
[KokkosParticleState.md].

## Example Input Syntax

!listing modules/dem/test/tests/kokkos/integration/ballistic.i block=UserObjects

!syntax parameters /UserObjects/KokkosParticleCloud

!syntax inputs /UserObjects/KokkosParticleCloud

!syntax children /UserObjects/KokkosParticleCloud
