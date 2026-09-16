#pragma once

#include "KokkosTypes.h"

#include "libmesh/point.h"
#include "libmesh/vector_value.h"

#include <vector>

namespace DEM
{

/**
 * Planar walls given analytically (plan layer L5), each by a point on it and its unit normal
 * pointing into the domain. The sphere-plane narrow phase is the signed distance of the particle
 * center to the plane, so a wall extends without bound and the mesh must keep the particle centers
 * on its inner side. Walls are rigid: the contact model sees an infinite-mass partner, so the
 * effective mass of a sphere-wall contact is the sphere's own. A wall may translate at a
 * velocity, prescribed or set by a servo, which its contacts see and which carries its point
 * along every substep.
 */
struct AnalyticWalls
{
  /// Number of walls
  std::size_t n = 0;
  /// A point on each wall
  ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight> points;
  /// Unit normal of each wall, pointing into the domain
  ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight> normals;
  /// Translation velocity of each wall
  ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight> velocities;
  /// Whether any wall moves
  bool moving = false;

  AnalyticWalls() = default;

  /// Build the walls from one point and one normal each; normals are normalized here
  AnalyticWalls(const std::vector<libMesh::Point> & wall_points,
                const std::vector<libMesh::RealVectorValue> & wall_normals)
    : n(wall_points.size()),
      points("dem_wall_points", n),
      normals("dem_wall_normals", n),
      velocities("dem_wall_velocities", n)
  {
    auto points_host = ::Kokkos::create_mirror_view(points);
    auto normals_host = ::Kokkos::create_mirror_view(normals);
    for (std::size_t w = 0; w < n; ++w)
    {
      const auto normal = wall_normals[w].unit();
      for (unsigned int c = 0; c < 3; ++c)
      {
        points_host(w, c) = wall_points[w](c);
        normals_host(w, c) = normal(c);
      }
    }
    ::Kokkos::deep_copy(points, points_host);
    ::Kokkos::deep_copy(normals, normals_host);
  }

  /// Unit normal of wall w
  KOKKOS_INLINE_FUNCTION Moose::Kokkos::Real3 normal(const std::size_t w) const
  {
    return Moose::Kokkos::Real3(normals(w, 0), normals(w, 1), normals(w, 2));
  }

  /// Velocity of wall w
  KOKKOS_INLINE_FUNCTION Moose::Kokkos::Real3 velocity(const std::size_t w) const
  {
    return Moose::Kokkos::Real3(velocities(w, 0), velocities(w, 1), velocities(w, 2));
  }

  /// Set the velocity of every wall (host)
  void setVelocities(const std::vector<libMesh::RealVectorValue> & v)
  {
    auto host = ::Kokkos::create_mirror_view(velocities);
    moving = false;
    for (std::size_t w = 0; w < n; ++w)
      for (unsigned int c = 0; c < 3; ++c)
      {
        host(w, c) = v[w](c);
        moving = moving || v[w](c) != 0;
      }
    ::Kokkos::deep_copy(velocities, host);
  }

  /// Carry the points along the velocities over a substep
  void advance(const Real dt) const
  {
    if (!moving)
      return;
    const auto p = points;
    const auto v = velocities;
    ::Kokkos::parallel_for(
        "dem_analytic_wall_advance", n, KOKKOS_LAMBDA(const std::size_t w) {
          for (unsigned int c = 0; c < 3; ++c)
            p(w, c) += dt * v(w, c);
        });
  }

  ///@{
  /// The points as a flat host vector, for checkpoints
  std::vector<Real> savePoints() const
  {
    auto host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, points);
    return std::vector<Real>(host.data(), host.data() + 3 * n);
  }
  void loadPoints(const std::vector<Real> & flat)
  {
    auto host = ::Kokkos::create_mirror_view(points);
    for (std::size_t w = 0; w < n; ++w)
      for (unsigned int c = 0; c < 3; ++c)
        host(w, c) = flat[3 * w + c];
    ::Kokkos::deep_copy(points, host);
  }
  ///@}

  /// Overlap of a sphere of radius r centered at x with wall w, positive in contact
  KOKKOS_INLINE_FUNCTION Real overlap(const std::size_t w,
                                      const Moose::Kokkos::Real3 & x,
                                      const Real r) const
  {
    const Moose::Kokkos::Real3 p(points(w, 0), points(w, 1), points(w, 2));
    return r - (x - p).dot_product(normal(w));
  }
};

} // namespace DEM
