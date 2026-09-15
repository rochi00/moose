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
 * History of a contact, stored once per pair in the orientation of the lower-ID particle: the
 * tangential spring displacement of the lower particle relative to the upper (of the particle
 * relative to the wall). The side of the upper particle negates it.
 */
struct PairState
{
  Real delta_t[3] = {0, 0, 0};
};

/// Device hash map of the contact histories (plan decision D6)
using PairStateMap = ::Kokkos::UnorderedMap<PairKey, PairState>;

} // namespace DEM
