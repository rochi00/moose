#pragma once

#include "KokkosMesh.h"

#include "libmesh/kokkos_storage_policy.h"
#include "libmesh/kokkos_storage.h"
#include "libmesh/kokkos_vector_ops.h"

namespace DEM
{

/**
 * Device-resident particle state (layer L0 of the DEM module plan).
 *
 * Vector quantities use the libMesh Kokkos storage policy so kernels can operate on entries
 * through libMesh::Kokkos::vector_ref / Point without a third vector type. LayoutRight keeps each
 * particle's components contiguous, which matches Moose::Kokkos::Real3 and LAMMPS-KOKKOS; the
 * policy is a single alias so it can be benchmarked against LayoutLeft later (plan Section 4.4).
 *
 * Element IDs are the Kokkos-MOOSE contiguous local element IDs (Moose::Kokkos::Mesh), not libMesh
 * IDs, and are therefore only meaningful on the rank that owns the particle.
 */
struct ParticleCloud
{
  using StoragePolicy = libMesh::Kokkos::layout_right_storage_policy;
  using VectorView = StoragePolicy::vector_view;

  /// Global particle ID, unique for the lifetime of the run (plan decision D2)
  ::Kokkos::View<int64_t *> gid;
  /// Position
  VectorView x;
  /// Velocity
  VectorView v;
  /// Contiguous local element containing the particle; libMesh::DofObject::invalid_id if unresolved
  ::Kokkos::View<ContiguousElementID *> elem;
  /// Rank that should receive the particle when it has walked into a ghost element;
  /// libMesh::DofObject::invalid_processor_id when the particle is local or lost
  ::Kokkos::View<processor_id_type *> target_rank;
  /// Number of face-neighbor hops the tracker took for this particle in the last step
  ::Kokkos::View<unsigned int *> hops;

  /// Number of live particles; entries [n, capacity()) are unused
  std::size_t n = 0;

  /// Number of live particles
  std::size_t size() const { return n; }
  /// Number of particles the views can hold
  std::size_t capacity() const { return gid.extent(0); }

  /// Allocate all views for a given capacity with no live particles (contents uninitialized)
  void allocate(const std::size_t capacity)
  {
    gid = ::Kokkos::View<int64_t *>("dem_gid", capacity);
    x = libMesh::Kokkos::make_vector_storage<StoragePolicy>("dem_x", capacity);
    v = libMesh::Kokkos::make_vector_storage<StoragePolicy>("dem_v", capacity);
    elem = ::Kokkos::View<ContiguousElementID *>("dem_elem", capacity);
    target_rank = ::Kokkos::View<processor_id_type *>("dem_target_rank", capacity);
    hops = ::Kokkos::View<unsigned int *>("dem_hops", capacity);
    n = 0;
  }

  /// Copy particle i of another cloud into slot j of this one
  KOKKOS_INLINE_FUNCTION void copyFrom(const ParticleCloud & other,
                                       const std::size_t i,
                                       const std::size_t j) const
  {
    gid(j) = other.gid(i);
    for (unsigned int c = 0; c < 3; ++c)
    {
      x(j, c) = other.x(i, c);
      v(j, c) = other.v(i, c);
    }
    elem(j) = other.elem(i);
    target_rank(j) = other.target_rank(i);
    hops(j) = other.hops(i);
  }
};

} // namespace DEM
