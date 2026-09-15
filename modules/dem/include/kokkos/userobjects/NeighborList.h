#pragma once

#include "ParticleCloud.h"

namespace DEM
{

/**
 * Uniform-grid broad phase producing a full CSR neighbor list (layer L2 of the DEM module plan):
 * particle i's neighbors are pairs[offsets(i) .. offsets(i + 1)), every pair appearing from both
 * sides so force kernels can sum a particle's contacts without atomics.
 *
 * Particles are binned on a grid of cell size 2 r_max + skin and counting-sorted by
 * bin * n + index, a unique key, so the neighbor order is deterministic for a given particle
 * order. Two particles are neighbors when their center distance is below r_i + r_j + skin. The
 * list stays valid until some particle has moved more than skin / 2 since the build.
 */
class NeighborList
{
public:
  /// Build the list for the live particles of a cloud within a bounding box
  void build(const ParticleCloud & cloud,
             const Moose::Kokkos::Real3 & lower,
             const Moose::Kokkos::Real3 & upper,
             const Real skin);

  /// Whether some particle has moved more than skin / 2 since the build, or the cloud size changed
  bool stale(const ParticleCloud & cloud) const;

  /// Number of neighbor pairs, each counted once
  std::size_t numPairs() const { return _num_pairs; }

  /// Neighbor offsets into pairs(), of size n + 1
  const ::Kokkos::View<std::size_t *> & offsets() const { return _offsets; }
  /// Concatenated neighbor indices
  const ::Kokkos::View<std::size_t *> & pairs() const { return _pairs; }

private:
  /// Number of particles the list was built for
  std::size_t _n = 0;
  /// Skin distance the list was built with
  Real _skin = 0;
  std::size_t _num_pairs = 0;
  ::Kokkos::View<std::size_t *> _offsets;
  ::Kokkos::View<std::size_t *> _pairs;
  /// Positions at the time of the build, for the staleness check
  ParticleCloud::VectorView _x_at_build;

  ///@{
  /// Work buffers, sized to the particle count at the last build
  ::Kokkos::View<uint64_t *> _keys;
  ::Kokkos::View<std::size_t *> _order;
  ::Kokkos::View<std::size_t *> _bin_offsets;
  ::Kokkos::View<std::size_t *> _counts;
  ///@}
};

} // namespace DEM
