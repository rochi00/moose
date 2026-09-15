#pragma once

#include "KokkosTypes.h"

#include <Kokkos_UnorderedMap.hpp>

namespace DEM
{

/**
 * Key of a contact's history (plan layer L3): the global IDs of the two particles in increasing
 * order, or the particle's global ID and, for a wall contact, the wall's index encoded as
 * -1 - wall. The same key is formed on every rank holding either particle, so the history of a
 * contact across a partition is kept in identical copies on both sides.
 */
struct PairKey
{
  int64_t lower;
  int64_t upper;

  KOKKOS_INLINE_FUNCTION static PairKey pair(const int64_t gid_i, const int64_t gid_j)
  {
    return gid_i < gid_j ? PairKey{gid_i, gid_j} : PairKey{gid_j, gid_i};
  }
  KOKKOS_INLINE_FUNCTION static PairKey wall(const int64_t gid, const std::size_t w)
  {
    return {gid, -1 - static_cast<int64_t>(w)};
  }
};

/**
 * History of a contact, stored once per pair: the tangential spring displacement, in the
 * orientation of the lower-ID particle (of the lower particle relative to the upper, or of the
 * particle relative to the wall; the side of the upper particle negates it), and the rolling
 * spring displacement, which is the same seen from either side.
 */
struct PairState
{
  Real delta_t[3] = {0, 0, 0};
  Real delta_r[3] = {0, 0, 0};

  /// Whether any spring is stretched
  KOKKOS_INLINE_FUNCTION bool stretched() const
  {
    for (unsigned int c = 0; c < 3; ++c)
      if (delta_t[c] != 0 || delta_r[c] != 0)
        return true;
    return false;
  }
};

/// Device hash map of the contact histories (plan decision D6)
using PairStateMap = ::Kokkos::UnorderedMap<PairKey, PairState>;

} // namespace DEM
