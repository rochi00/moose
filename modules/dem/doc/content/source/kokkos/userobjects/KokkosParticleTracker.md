# KokkosParticleTracker

!syntax description /UserObjects/KokkosParticleTracker

## Overview

`KokkosParticleTracker` holds a cloud of particles on device and, every time it executes, advances
each particle by `velocity * dt` and re-resolves the local element containing it. Particle state is
stored in Kokkos views using the libMesh Kokkos storage policy (`libmesh/kokkos_storage_policy.h`).

Element tracking is incremental: each particle remembers its element, and after moving it is tested
against the outward face planes cached by the Kokkos mesh (`Moose::Kokkos::Mesh::getSideNormal`,
`getSideCentroid`). A particle in front of a face plane leaves through the face it is furthest in
front of and the test repeats on the face neighbor, up to `max_hops` times. This is exact only for
elements with planar faces.

The walk ends in one of three states:

- inside a local element;
- inside a one-layer ghost element, which has no cached face geometry. The particle is then sent
  to the rank owning that element, which resumes the walk from there. Rounds of walking and
  exchanging repeat until no particle is left in a ghost element, so a particle may cross several
  partitions in one step. Departing particles are found, compacted out of the cloud, and packed on
  device, so only the departing and arriving particles cross to the host for the TIMPI exchange;
  GPU-aware MPI is deferred;
- lost, because it walked out of the mesh or exhausted `max_hops`.

With `verify = true`, each step the device assignment is checked on host against libMesh's
`PointLocator` and any mismatch is an error.

Counters are exposed through [KokkosParticleTrackerValue.md].

## Example Input Syntax

!listing test/tests/kokkos/tracker/walk_2d.i block=UserObjects

!syntax parameters /UserObjects/KokkosParticleTracker

!syntax inputs /UserObjects/KokkosParticleTracker

!syntax children /UserObjects/KokkosParticleTracker
