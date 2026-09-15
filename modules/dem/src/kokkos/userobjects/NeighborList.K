// Kokkos-MOOSE's KokkosHeader.h defines macros (MemSpace, KOKKOS_IF_ON_HOST) that break the Kokkos
// sort headers, so they must come first and are only included in this translation unit
#include <Kokkos_Sort.hpp>
#ifdef DEM_HAVE_ARBORX
#include <ArborX_LinearBVH.hpp>
#endif

#include "NeighborList.h"
#include "MooseError.h"

#ifdef DEM_HAVE_ARBORX
namespace DEM
{
/// The first n spheres of a cloud as ArborX primitives (their bounding boxes) or, inflated by the
/// skin, as intersection predicates carrying the particle index
struct SphereBoxes
{
  ParticleCloud cloud;
  std::size_t n;
  Real inflation;
};

/// Keeps a candidate pair when the exact criterion holds, distance below r_i + r_j + skin
struct SphereFilter
{
  ParticleCloud cloud;
  Real skin;

  template <typename Predicate, typename Value, typename Output>
  KOKKOS_FUNCTION void
  operator()(const Predicate & predicate, const Value & value, const Output & out) const
  {
    const std::size_t i = ArborX::getData(predicate);
    const std::size_t j = value.index;
    if (j == i)
      return;
    const Moose::Kokkos::Real3 d(
        cloud.x(i, 0) - cloud.x(j, 0), cloud.x(i, 1) - cloud.x(j, 1), cloud.x(i, 2) - cloud.x(j, 2));
    if (d.norm() < cloud.r(i) + cloud.r(j) + skin)
      out(static_cast<int>(j));
  }
};
} // namespace DEM

template <>
struct ArborX::AccessTraits<DEM::SphereBoxes>
{
  using memory_space = ::Kokkos::DefaultExecutionSpace::memory_space;
  static KOKKOS_FUNCTION std::size_t size(const DEM::SphereBoxes & s) { return s.n; }
  static KOKKOS_FUNCTION auto get(const DEM::SphereBoxes & s, const std::size_t i)
  {
    const Real h = s.cloud.r(i) + s.inflation;
    return ArborX::Box<3, Real>{{s.cloud.x(i, 0) - h, s.cloud.x(i, 1) - h, s.cloud.x(i, 2) - h},
                                {s.cloud.x(i, 0) + h, s.cloud.x(i, 1) + h, s.cloud.x(i, 2) + h}};
  }
};
#endif

namespace DEM
{

void
NeighborList::build(const ParticleCloud & cloud,
                    const std::size_t n_local,
                    const std::size_t n_total,
                    const Moose::Kokkos::Real3 & lower,
                    const Moose::Kokkos::Real3 & upper,
                    const Real skin,
                    const BroadPhase broad_phase)
{
  const auto n = n_total;
  _n_local = n_local;
  _skin = skin;

  // A rank without particles still needs a valid (empty) list
  if (n == 0)
  {
    if (_offsets.extent(0) < 1)
      _offsets = ::Kokkos::View<std::size_t *>("dem_nl_offsets", 1);
    ::Kokkos::deep_copy(_offsets, 0);
    _num_pairs = 0;
    _num_interior = _num_boundary = 0;
    return;
  }

  if (_offsets.extent(0) < n + 1)
    _offsets = ::Kokkos::View<std::size_t *>("dem_nl_offsets", n + 1);
  if (_x_at_build.extent(0) < n_local)
    _x_at_build = libMesh::Kokkos::make_vector_storage<ParticleCloud::StoragePolicy>(
        "dem_nl_x_at_build", n_local);
  const auto x_at_build = _x_at_build;
  ::Kokkos::parallel_for(
      "dem_nl_x_at_build", n_local, KOKKOS_LAMBDA(const std::size_t i) {
        for (unsigned int c = 0; c < 3; ++c)
          x_at_build(i, c) = cloud.x(i, c);
      });

  if (broad_phase == BroadPhase::ARBORX_BVH)
    buildArborX(cloud, n, skin);
  else
    buildUniformGrid(cloud, n, lower, upper, skin);
  sortNeighbors(n);
  const auto offsets = _offsets;
  const auto pairs = _pairs;

  // Every pair is listed from both sides: count a local pair from its lower index and a pair with
  // a ghost from the rank owning the lower global ID
  std::size_t num_pairs = 0;
  ::Kokkos::parallel_reduce(
      "dem_nl_num_pairs",
      n_local,
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum) {
        for (auto k = offsets(i); k < offsets(i + 1); ++k)
        {
          const auto j = pairs(k);
          sum += j < n_local ? i < j : cloud.gid(i) < cloud.gid(j);
        }
      },
      num_pairs);
  _num_pairs = num_pairs;

  // Split the local particles by whether they have a ghost neighbor, compacting each side with a
  // scan; the boundary list is the interior list's complement
  if (_interior.extent(0) < n_local)
  {
    _interior = ::Kokkos::View<std::size_t *>("dem_nl_interior", n_local);
    _boundary = ::Kokkos::View<std::size_t *>("dem_nl_boundary", n_local);
  }
  const auto interior = _interior;
  const auto boundary = _boundary;
  std::size_t num_interior = 0;
  ::Kokkos::parallel_scan(
      "dem_nl_split",
      n_local,
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & n_interior_before, const bool final) {
        bool is_interior = true;
        for (auto k = offsets(i); k < offsets(i + 1) && is_interior; ++k)
          is_interior = pairs(k) < n_local;
        if (final)
        {
          if (is_interior)
            interior(n_interior_before) = i;
          else
            boundary(i - n_interior_before) = i;
        }
        n_interior_before += is_interior;
      },
      num_interior);
  _num_interior = num_interior;
  _num_boundary = n_local - num_interior;
}

void
NeighborList::buildUniformGrid(const ParticleCloud & cloud,
                               const std::size_t n,
                               const Moose::Kokkos::Real3 & lower,
                               const Moose::Kokkos::Real3 & upper,
                               const Real skin)
{
  // Copies for capture by value in the kernels
  const Moose::Kokkos::Real3 lo = lower;

  // Cell size from the largest radius so any two overlapping-plus-skin spheres are in the same or
  // adjacent cells
  Real r_max = 0;
  ::Kokkos::parallel_reduce(
      "dem_nl_rmax",
      n,
      KOKKOS_LAMBDA(const std::size_t i, Real & max) { max = cloud.r(i) > max ? cloud.r(i) : max; },
      ::Kokkos::Max<Real>(r_max));
  const Real h = 2 * r_max + skin;
  int dims[3];
  for (const auto c : make_range(3))
    dims[c] = std::max(1, static_cast<int>(std::floor((upper(c) - lower(c)) / h)));
  const int nx = dims[0], ny = dims[1], nz = dims[2];
  const std::size_t num_bins = static_cast<std::size_t>(nx) * ny * nz;

  if (_keys.extent(0) < n)
  {
    _keys = ::Kokkos::View<uint64_t *>("dem_nl_keys", n);
    _order = ::Kokkos::View<std::size_t *>("dem_nl_order", n);
  }
  if (_counts.extent(0) < num_bins)
    _counts = ::Kokkos::View<std::size_t *>("dem_nl_counts", num_bins);
  if (_bin_offsets.extent(0) < num_bins + 1)
    _bin_offsets = ::Kokkos::View<std::size_t *>("dem_nl_bin_offsets", num_bins + 1);
  const auto keys = _keys;
  const auto order = _order;
  const auto offsets = _offsets;
  const auto counts = _counts;
  const auto bin_offsets = _bin_offsets;

  // Bin index of each particle, clamped so particles outside the box land in edge bins
  auto bin_of = KOKKOS_LAMBDA(const std::size_t i)
  {
    int b[3];
    for (unsigned int c = 0; c < 3; ++c)
    {
      const int idx = static_cast<int>(::Kokkos::floor((cloud.x(i, c) - lo(c)) / h));
      b[c] = idx < 0 ? 0 : (idx >= dims[c] ? dims[c] - 1 : idx);
    }
    return static_cast<std::size_t>(b[0]) + nx * (static_cast<std::size_t>(b[1]) + ny * b[2]);
  };

  // Counting sort by (bin, index): unique keys make the order deterministic
  ::Kokkos::deep_copy(::Kokkos::subview(counts, std::make_pair(std::size_t(0), num_bins)), 0);
  ::Kokkos::parallel_for(
      "dem_nl_keys", n, KOKKOS_LAMBDA(const std::size_t i) {
        const auto bin = bin_of(i);
        keys(i) = static_cast<uint64_t>(bin) * n + i;
        order(i) = i;
        ::Kokkos::atomic_inc(&counts(bin));
      });
  ::Kokkos::Experimental::sort_by_key(::Kokkos::DefaultExecutionSpace(),
                                      ::Kokkos::subview(keys, std::make_pair(std::size_t(0), n)),
                                      ::Kokkos::subview(order, std::make_pair(std::size_t(0), n)));
  ::Kokkos::parallel_scan(
      "dem_nl_bin_scan",
      num_bins + 1,
      KOKKOS_LAMBDA(const std::size_t b, std::size_t & sum, const bool final) {
        if (final)
          bin_offsets(b) = sum;
        if (b < num_bins)
          sum += counts(b);
      });

  // Neighbors of particle i: every particle in the 27 surrounding bins closer than r_i + r_j + skin
  auto for_each_neighbor = KOKKOS_LAMBDA(const std::size_t i, auto && action)
  {
    const Moose::Kokkos::Real3 xi(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2));
    int b[3];
    for (unsigned int c = 0; c < 3; ++c)
    {
      const int idx = static_cast<int>(::Kokkos::floor((cloud.x(i, c) - lo(c)) / h));
      b[c] = idx < 0 ? 0 : (idx >= dims[c] ? dims[c] - 1 : idx);
    }
    for (int dz = -1; dz <= 1; ++dz)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx)
        {
          const int bx = b[0] + dx, by = b[1] + dy, bz = b[2] + dz;
          if (bx < 0 || bx >= nx || by < 0 || by >= ny || bz < 0 || bz >= nz)
            continue;
          const std::size_t bin =
              static_cast<std::size_t>(bx) + nx * (static_cast<std::size_t>(by) + ny * bz);
          for (auto k = bin_offsets(bin); k < bin_offsets(bin + 1); ++k)
          {
            const auto j = order(k);
            if (j == i)
              continue;
            const Moose::Kokkos::Real3 xj(cloud.x(j, 0), cloud.x(j, 1), cloud.x(j, 2));
            const Real cutoff = cloud.r(i) + cloud.r(j) + skin;
            if ((xi - xj).norm() < cutoff)
              action(j);
          }
        }
  };

  // Count, scan, fill
  ::Kokkos::parallel_for(
      "dem_nl_count", n, KOKKOS_LAMBDA(const std::size_t i) {
        std::size_t count = 0;
        for_each_neighbor(i, [&](const std::size_t) { ++count; });
        offsets(i) = count;
      });
  std::size_t total = 0;
  ::Kokkos::parallel_scan(
      "dem_nl_offset_scan",
      n + 1,
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum, const bool final) {
        const auto count = i < n ? offsets(i) : 0;
        if (final)
          offsets(i) = sum;
        sum += count;
      },
      total);
  if (_pairs.extent(0) < total)
    _pairs = ::Kokkos::View<std::size_t *>("dem_nl_pairs", total);
  const auto pairs = _pairs;
  ::Kokkos::parallel_for(
      "dem_nl_fill", n, KOKKOS_LAMBDA(const std::size_t i) {
        auto k = offsets(i);
        for_each_neighbor(i, [&](const std::size_t j) { pairs(k++) = j; });
      });
}

void
NeighborList::buildArborX(const ParticleCloud & cloud, const std::size_t n, const Real skin)
{
#ifdef DEM_HAVE_ARBORX
  // A hierarchy over the spheres' bounding boxes, queried with each sphere's box inflated by the
  // skin; box-box overlap is a superset of the sphere criterion, which the callback then applies
  // exactly. ArborX returns the CRS list in its own integer type, copied into ours
  const ::Kokkos::DefaultExecutionSpace space;
  ArborX::BoundingVolumeHierarchy bvh(
      space, ArborX::Experimental::attach_indices<int>(SphereBoxes{cloud, n, 0}));
  bvh.query(space,
            ArborX::Experimental::attach_indices<int>(
                ArborX::Experimental::make_intersects(SphereBoxes{cloud, n, skin})),
            SphereFilter{cloud, skin},
            _arborx_pairs,
            _arborx_offsets);

  const auto arborx_offsets = _arborx_offsets;
  const auto arborx_pairs = _arborx_pairs;
  const std::size_t total = arborx_pairs.extent(0);
  if (_pairs.extent(0) < total)
    _pairs = ::Kokkos::View<std::size_t *>("dem_nl_pairs", total);
  const auto offsets = _offsets;
  const auto pairs = _pairs;
  ::Kokkos::parallel_for(
      "dem_nl_arborx_offsets", n + 1, KOKKOS_LAMBDA(const std::size_t i) {
        offsets(i) = arborx_offsets(i);
      });
  ::Kokkos::parallel_for(
      "dem_nl_arborx_pairs", total, KOKKOS_LAMBDA(const std::size_t k) {
        pairs(k) = arborx_pairs(k);
      });
#else
  libmesh_ignore(cloud, n, skin);
  mooseError("The ArborX broad phase needs the ArborX submodule (modules/dem/contrib/arborx)");
#endif
}

void
NeighborList::sortNeighbors(const std::size_t n)
{
  // Insertion sort of each particle's few neighbors, in place
  const auto offsets = _offsets;
  const auto pairs = _pairs;
  ::Kokkos::parallel_for(
      "dem_nl_sort", n, KOKKOS_LAMBDA(const std::size_t i) {
        for (auto k = offsets(i) + 1; k < offsets(i + 1); ++k)
        {
          const auto j = pairs(k);
          auto m = k;
          for (; m > offsets(i) && pairs(m - 1) > j; --m)
            pairs(m) = pairs(m - 1);
          pairs(m) = j;
        }
      });
}

bool
NeighborList::stale(const ParticleCloud & cloud, const std::size_t n_local) const
{
  if (n_local != _n_local)
    return true;

  const auto x_at_build = _x_at_build;
  Real max_displacement_sq = 0;
  ::Kokkos::parallel_reduce(
      "dem_nl_stale",
      n_local,
      KOKKOS_LAMBDA(const std::size_t i, Real & max) {
        Real d_sq = 0;
        for (unsigned int c = 0; c < 3; ++c)
          d_sq += (cloud.x(i, c) - x_at_build(i, c)) * (cloud.x(i, c) - x_at_build(i, c));
        max = d_sq > max ? d_sq : max;
      },
      ::Kokkos::Max<Real>(max_displacement_sq));
  // Two particles approaching each other can each cover skin / 2 before an unlisted pair touches
  return max_displacement_sq > 0.25 * _skin * _skin;
}

} // namespace DEM
