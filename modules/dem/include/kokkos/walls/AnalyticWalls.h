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
  /// Rigid rotation of the whole set of walls about an axis through a center at a constant
  /// angular velocity (rad per unit time; zero when not rotating)
  Moose::Kokkos::Real3 rotation_center = Moose::Kokkos::Real3(0);
  Moose::Kokkos::Real3 angular_velocity = Moose::Kokkos::Real3(0);
  bool rotating = false;

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

  /// Velocity of wall w at the point x of its surface: the translation plus the rotation's
  /// omega x (x - center)
  KOKKOS_INLINE_FUNCTION Moose::Kokkos::Real3 velocity(const std::size_t w,
                                                       const Moose::Kokkos::Real3 & x) const
  {
    Moose::Kokkos::Real3 v(velocities(w, 0), velocities(w, 1), velocities(w, 2));
    if (rotating)
      v += angular_velocity.cross_product(x - rotation_center);
    return v;
  }

  /// Set the rigid rotation of the walls (host)
  void setRotation(const libMesh::Point & center, const libMesh::RealVectorValue & omega)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      rotation_center(c) = center(c);
      angular_velocity(c) = omega(c);
    }
    rotating = omega.norm() > 0;
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

  /// Carry the points along the velocities over a substep, and turn the points and normals by
  /// the rotation over it (exactly, by Rodrigues' formula)
  void advance(const Real dt) const
  {
    if (!moving && !rotating)
      return;
    const auto p = points;
    const auto nrm = normals;
    const auto v = velocities;
    const auto center = rotation_center;
    const auto omega = angular_velocity;
    const bool rotate = rotating;
    ::Kokkos::parallel_for(
        "dem_analytic_wall_advance", n, KOKKOS_LAMBDA(const std::size_t w) {
          for (unsigned int c = 0; c < 3; ++c)
            p(w, c) += dt * v(w, c);
          if (!rotate)
            return;
          const Real theta = omega.norm() * dt;
          const Moose::Kokkos::Real3 axis = (1.0 / omega.norm()) * omega;
          const Real cos_t = ::Kokkos::cos(theta), sin_t = ::Kokkos::sin(theta);
          Moose::Kokkos::Real3 r(p(w, 0) - center(0), p(w, 1) - center(1), p(w, 2) - center(2));
          Moose::Kokkos::Real3 m(nrm(w, 0), nrm(w, 1), nrm(w, 2));
          r = cos_t * r + sin_t * axis.cross_product(r) + (1 - cos_t) * axis.dot_product(r) * axis;
          m = cos_t * m + sin_t * axis.cross_product(m) + (1 - cos_t) * axis.dot_product(m) * axis;
          for (unsigned int c = 0; c < 3; ++c)
          {
            p(w, c) = center(c) + r(c);
            nrm(w, c) = m(c);
          }
        });
  }

  ///@{
  /// The points and normals as one flat host vector, for checkpoints
  std::vector<Real> savePoints() const
  {
    auto host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, points);
    auto host_n = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, normals);
    std::vector<Real> flat(host.data(), host.data() + 3 * n);
    flat.insert(flat.end(), host_n.data(), host_n.data() + 3 * n);
    return flat;
  }
  void loadPoints(const std::vector<Real> & flat)
  {
    auto host = ::Kokkos::create_mirror_view(points);
    auto host_n = ::Kokkos::create_mirror_view(normals);
    for (std::size_t w = 0; w < n; ++w)
      for (unsigned int c = 0; c < 3; ++c)
      {
        host(w, c) = flat[3 * w + c];
        host_n(w, c) = flat[3 * n + 3 * w + c];
      }
    ::Kokkos::deep_copy(points, host);
    ::Kokkos::deep_copy(normals, host_n);
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
