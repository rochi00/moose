#pragma once

#include "KokkosTypes.h"
#include "KokkosMesh.h"

#include "libmesh/bounding_box.h"

#include <vector>
#include <array>

class MooseMesh;

namespace DEM
{

/// A sphere's contact with a wall face: the closest point of the face, the unit normal from it
/// toward the sphere center, the overlap, the face, and the index of the face's boundary, which
/// keys the contact history
struct WallContact
{
  Moose::Kokkos::Real3 point;
  Moose::Kokkos::Real3 normal;
  Real overlap;
  std::size_t face;
  std::size_t boundary;
  /// Velocity of the wall at the point (zero unless the walls move)
  Moose::Kokkos::Real3 velocity;
};

/**
 * Walls derived from mesh sidesets (plan layer L5, decision D9): the faces of the selected
 * boundaries of the local and one-layer ghost elements, triangulated (a quad split along its
 * first diagonal; in 2D a face is an edge segment), exported to device with their inward
 * normals and a global ID. Each element (local or ghost) lists the faces whose bounding box,
 * inflated by the reach, meets its own, so a particle only tests the faces near the element it
 * is tracked in.
 *
 * The sphere-face narrow phase finds the closest point of each candidate face (Ericson 2005).
 * Contacts are then reduced so that a sphere over an edge or vertex shared by several faces
 * gets one force: contacts with the same closest point are merged (lowest face kept) and an
 * edge or vertex contact is dropped when its point lies on a face that has an interior contact,
 * so a triangulated plane acts exactly as a plane and a convex edge or vertex as a single
 * contact along the line to the center.
 *
 * The contact history is keyed by the face's boundary, not the face, so it carries across the
 * faces of one sideset as the contact point slides over them, a triangulated plane keeping its
 * tangential spring like a plane. Faces a particle can touch at the same time with different
 * normals (the floor and a side of a box) should therefore be in different sidesets.
 *
 * The walls can move: given new vertex positions once per MOOSE step, each vertex is assigned
 * the velocity that carries it there over the step, and the faces advance with it every substep,
 * the velocity of a contact point being interpolated from the vertices. The candidate lists are
 * not rebuilt, so the motion must stay within the reach the lists were built with.
 */
struct SidesetWalls
{
  /// Number of faces held (local and ghost elements' sideset faces)
  std::size_t n = 0;
  /// Vertices of each face, three per face (the third repeats the second for a segment)
  ::Kokkos::View<Real * [3][3], ::Kokkos::LayoutRight> vertices;
  /// Number of vertices of each face: 3 for a triangle, 2 for a segment
  ::Kokkos::View<unsigned int *> num_vertices;
  /// Unit normal of each face pointing into the domain
  ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight> normals;
  /// Velocity of each vertex of each face, zero unless the walls move
  ::Kokkos::View<Real * [3][3], ::Kokkos::LayoutRight> velocities;
  /// Mesh node of each vertex of each face, on host, to sample the displacement at
  std::vector<std::array<libMesh::dof_id_type, 3>> nodes;
  /// Whether any vertex velocity is nonzero
  bool moving = false;
  /// Index of each face's boundary in the list the walls were built from
  ::Kokkos::View<std::size_t *> boundary;
  /// Number of boundaries the walls were built from
  std::size_t num_boundaries = 0;
  ///@{
  /// CSR list of the candidate faces of each local and ghost element, by contiguous element ID
  ::Kokkos::View<std::size_t *> elem_offsets;
  ::Kokkos::View<std::size_t *> elem_faces;
  ///@}
  /// Most contacts a sphere is evaluated against at once
  static constexpr unsigned int max_contacts = 16;
  /// Set by contacts() when a sphere had more candidate contacts than max_contacts; checked
  /// after the kernels so the truncation is an error rather than a silent loss of contacts
  ::Kokkos::View<int> overflow;

  /**
   * Build the walls from the given boundaries of the local and one-layer ghost elements
   * @param reach The distance within which a face is a candidate of an element (largest radius
   *        plus the neighbor-list skin)
   */
  void build(const MooseMesh & mesh,
             const Moose::Kokkos::Mesh & kokkos_mesh,
             const std::vector<const libMesh::Elem *> & cid_to_elem,
             const std::vector<BoundaryID> & boundaries,
             const Real reach);

  /**
   * Set the vertex velocities so that the vertices reach the given positions (three per face,
   * in face order) over a time dt, and recompute the normals for the positions reached; with a
   * zero dt the vertices are placed there at once, at rest
   */
  void move(const std::vector<Point> & positions, const Real dt);
  /// Advance the vertices by their velocities over a substep
  void advance(const Real dt);
  ///@{
  /// Copy the vertex positions and velocities out to, and back in from, flat host vectors
  /// (nine per face each), for checkpoints
  void save(std::vector<Real> & positions, std::vector<Real> & velocities) const;
  void load(const std::vector<Real> & positions, const std::vector<Real> & velocities);
  ///@}

  /// Closest point of face f to x, and whether it is in the face's interior (not on an edge or
  /// vertex)
  KOKKOS_INLINE_FUNCTION Moose::Kokkos::Real3
  closestPoint(const std::size_t f, const Moose::Kokkos::Real3 & x, bool & interior) const;
  /// Velocity of face f at a point on it, interpolated from its vertices
  KOKKOS_INLINE_FUNCTION Moose::Kokkos::Real3 velocityAt(const std::size_t f,
                                                         const Moose::Kokkos::Real3 & p) const;

  /**
   * The contacts of a sphere in element elem with the walls, reduced as described above
   * @param reach Faces closer than this to the center are contacts (the radius for forces, the
   *        radius plus the skin to list the histories to keep)
   * @param out Filled with at most max_contacts contacts; the overlap is r minus the distance
   *        along the line from the closest point, or the signed distance to the face plane for
   *        an interior contact, so a center behind the face is pushed back through it
   * @returns The number of contacts
   */
  KOKKOS_INLINE_FUNCTION unsigned int contacts(const ContiguousElementID elem,
                                               const Moose::Kokkos::Real3 & x,
                                               const Real r,
                                               const Real reach,
                                               WallContact * const out) const;
};

KOKKOS_INLINE_FUNCTION Moose::Kokkos::Real3
SidesetWalls::closestPoint(const std::size_t f, const Moose::Kokkos::Real3 & x, bool & interior) const
{
  const Moose::Kokkos::Real3 a(vertices(f, 0, 0), vertices(f, 0, 1), vertices(f, 0, 2));
  const Moose::Kokkos::Real3 b(vertices(f, 1, 0), vertices(f, 1, 1), vertices(f, 1, 2));
  interior = false;

  if (num_vertices(f) == 2)
  {
    // Segment: clamp the projection onto it
    const Moose::Kokkos::Real3 ab = b - a;
    const Real t = (x - a).dot_product(ab) / ab.dot_product(ab);
    if (t <= 0)
      return a;
    if (t >= 1)
      return b;
    interior = true;
    return a + t * ab;
  }

  // Triangle, Ericson's "Real-Time Collision Detection" 5.1.5: test the Voronoi regions of the
  // vertices, then the edges, else the face
  const Moose::Kokkos::Real3 c(vertices(f, 2, 0), vertices(f, 2, 1), vertices(f, 2, 2));
  const Moose::Kokkos::Real3 ab = b - a, ac = c - a, ax = x - a;
  const Real d1 = ab.dot_product(ax), d2 = ac.dot_product(ax);
  if (d1 <= 0 && d2 <= 0)
    return a;
  const Moose::Kokkos::Real3 bx = x - b;
  const Real d3 = ab.dot_product(bx), d4 = ac.dot_product(bx);
  if (d3 >= 0 && d4 <= d3)
    return b;
  const Real vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0)
    return a + (d1 / (d1 - d3)) * ab;
  const Moose::Kokkos::Real3 cx = x - c;
  const Real d5 = ab.dot_product(cx), d6 = ac.dot_product(cx);
  if (d6 >= 0 && d5 <= d6)
    return c;
  const Real vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0)
    return a + (d2 / (d2 - d6)) * ac;
  const Real va = d3 * d6 - d5 * d4;
  if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0)
    return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b);
  const Real denom = 1.0 / (va + vb + vc);
  interior = true;
  return a + (vb * denom) * ab + (vc * denom) * ac;
}

KOKKOS_INLINE_FUNCTION Moose::Kokkos::Real3
SidesetWalls::velocityAt(const std::size_t f, const Moose::Kokkos::Real3 & p) const
{
  if (!moving)
    return Moose::Kokkos::Real3(0);
  const Moose::Kokkos::Real3 a(vertices(f, 0, 0), vertices(f, 0, 1), vertices(f, 0, 2));
  const Moose::Kokkos::Real3 b(vertices(f, 1, 0), vertices(f, 1, 1), vertices(f, 1, 2));
  const Moose::Kokkos::Real3 va(velocities(f, 0, 0), velocities(f, 0, 1), velocities(f, 0, 2));
  const Moose::Kokkos::Real3 vb(velocities(f, 1, 0), velocities(f, 1, 1), velocities(f, 1, 2));
  if (num_vertices(f) == 2)
  {
    const Moose::Kokkos::Real3 ab = b - a;
    const Real t = (p - a).dot_product(ab) / ab.dot_product(ab);
    return (1 - t) * va + t * vb;
  }
  // Barycentric weights of p (Ericson 3.4)
  const Moose::Kokkos::Real3 c(vertices(f, 2, 0), vertices(f, 2, 1), vertices(f, 2, 2));
  const Moose::Kokkos::Real3 vc(velocities(f, 2, 0), velocities(f, 2, 1), velocities(f, 2, 2));
  const Moose::Kokkos::Real3 v0 = b - a, v1 = c - a, v2 = p - a;
  const Real d00 = v0.dot_product(v0), d01 = v0.dot_product(v1), d11 = v1.dot_product(v1),
             d20 = v2.dot_product(v0), d21 = v2.dot_product(v1);
  const Real denom = d00 * d11 - d01 * d01;
  const Real wb = (d11 * d20 - d01 * d21) / denom, wc = (d00 * d21 - d01 * d20) / denom;
  return (1 - wb - wc) * va + wb * vb + wc * vc;
}

KOKKOS_INLINE_FUNCTION unsigned int
SidesetWalls::contacts(const ContiguousElementID elem,
                       const Moose::Kokkos::Real3 & x,
                       const Real r,
                       const Real reach,
                       WallContact * const out) const
{
  if (n == 0 || elem >= elem_offsets.extent(0) - 1)
    return 0;
  unsigned int count = 0;
  bool interior_flags[max_contacts];
  std::size_t faces[max_contacts];
  for (auto k = elem_offsets(elem); k < elem_offsets(elem + 1); ++k)
  {
    const auto f = elem_faces(k);
    bool interior;
    const auto p = closestPoint(f, x, interior);
    const Moose::Kokkos::Real3 inward(normals(f, 0), normals(f, 1), normals(f, 2));
    const Moose::Kokkos::Real3 d = x - p;
    Real overlap;
    Moose::Kokkos::Real3 normal;
    if (interior)
    {
      // Signed: a center behind the face is still pushed back through it
      const Real height = d.dot_product(inward);
      overlap = r - height;
      normal = inward;
      if (height > reach)
        continue;
    }
    else
    {
      const Real distance = d.norm();
      if (distance >= reach || distance == 0)
        continue;
      overlap = r - distance;
      normal = (1.0 / distance) * d;
    }
    if (count == max_contacts)
    {
      overflow() = 1;
      break;
    }
    out[count] = {p, normal, overlap, f, boundary(f), velocityAt(f, p)};
    interior_flags[count] = interior;
    faces[count] = f;
    ++count;
  }

  // Reduce: merge contacts at the same point (keep the lowest ID), and drop an edge or vertex
  // contact whose point lies on a face with an interior contact
  const Real tol = 1e-12 * (r > 0 ? r : 1);
  bool drop[max_contacts];
  for (unsigned int i = 0; i < count; ++i)
  {
    drop[i] = false;
    for (unsigned int j = 0; j < count && !drop[i]; ++j)
    {
      if (j == i)
        continue;
      if ((out[i].point - out[j].point).norm() <= tol)
        drop[i] = out[j].face < out[i].face;
      else if (!interior_flags[i] && interior_flags[j])
      {
        bool unused;
        drop[i] = (closestPoint(faces[j], out[i].point, unused) - out[i].point).norm() <= tol;
      }
    }
  }
  unsigned int kept = 0;
  for (unsigned int i = 0; i < count; ++i)
    if (!drop[i])
      out[kept++] = out[i];
  return kept;
}

} // namespace DEM
