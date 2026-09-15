#include "SidesetWalls.h"

#include "MooseMesh.h"

#include "libmesh/elem.h"

#include <map>
#include <algorithm>

namespace DEM
{

void
SidesetWalls::build(const MooseMesh & mesh,
                    const Moose::Kokkos::Mesh & kokkos_mesh,
                    const std::vector<const libMesh::Elem *> & cid_to_elem,
                    const std::vector<BoundaryID> & boundaries,
                    const Real reach)
{
  libmesh_ignore(kokkos_mesh);
  const auto num_elems = cid_to_elem.size();
  num_boundaries = boundaries.size();

  // Faces on host: vertices, vertex count, inward normal, boundary index, bounding box
  struct Face
  {
    Point v[3];
    unsigned int nv;
    Point normal;
    std::size_t boundary;
    libMesh::BoundingBox box;
  };
  std::vector<Face> faces;
  const auto add = [&faces](const Point & a,
                            const Point & b,
                            const Point & c,
                            const unsigned int nv,
                            const Point & centroid,
                            const std::size_t boundary)
  {
    Face face{{a, b, c}, nv, Point(), boundary, libMesh::BoundingBox()};
    // Normal toward the element's centroid, which is inside the domain
    face.normal = nv == 3 ? (b - a).cross(c - a) : Point(-(b - a)(1), (b - a)(0), 0);
    face.normal /= face.normal.norm();
    if (face.normal * (centroid - a) < 0)
      face.normal *= -1;
    for (const auto k : make_range(nv))
      face.box.union_with(face.v[k]);
    faces.push_back(face);
  };
  for (const auto cid : make_range(num_elems))
  {
    const auto * elem = cid_to_elem[cid];
    for (const auto s : make_range(elem->n_sides()))
    {
      // A side on several of the boundaries goes with the first listed
      const auto ids = mesh.getBoundaryIDs(elem, s);
      std::size_t b = 0;
      for (; b < boundaries.size(); ++b)
        if (std::find(ids.begin(), ids.end(), boundaries[b]) != ids.end())
          break;
      if (b == boundaries.size())
        continue;
      const auto side = elem->side_ptr(s);
      const auto centroid = elem->vertex_average();
      // Vertices only: higher-order sides are used through their corner nodes
      std::vector<Point> v;
      for (const auto k : make_range(side->n_vertices()))
        v.push_back(side->point(k));
      if (v.size() == 2)
        add(v[0], v[1], v[1], 2, centroid, b);
      else if (v.size() == 3)
        add(v[0], v[1], v[2], 3, centroid, b);
      else if (v.size() == 4)
      {
        add(v[0], v[1], v[2], 3, centroid, b);
        add(v[0], v[2], v[3], 3, centroid, b);
      }
      else
        mooseError("Unsupported side with ", v.size(), " vertices for a sideset wall");
    }
  }
  n = faces.size();

  // Candidate faces of each element: those whose bounding box, inflated by the reach, meets
  // the element's. Faces are binned on a host grid so the search is linear in their number
  Real cell = 2 * reach;
  for (const auto & face : faces)
    for (const auto c : make_range(3))
      cell = std::max(cell, face.box.max()(c) - face.box.min()(c) + 2 * reach);
  std::map<std::array<long, 3>, std::vector<std::size_t>> bins;
  const auto cell_of = [cell](const Point & p)
  {
    return std::array<long, 3>{static_cast<long>(std::floor(p(0) / cell)),
                               static_cast<long>(std::floor(p(1) / cell)),
                               static_cast<long>(std::floor(p(2) / cell))};
  };
  for (const auto f : index_range(faces))
  {
    const auto lo = cell_of(faces[f].box.min()), hi = cell_of(faces[f].box.max());
    for (long i = lo[0]; i <= hi[0]; ++i)
      for (long j = lo[1]; j <= hi[1]; ++j)
        for (long k = lo[2]; k <= hi[2]; ++k)
          bins[{i, j, k}].push_back(f);
  }
  std::vector<std::size_t> offsets_host(num_elems + 1, 0), faces_host;
  for (const auto cid : make_range(num_elems))
  {
    auto box = cid_to_elem[cid]->loose_bounding_box();
    box.first -= Point(reach, reach, reach);
    box.second += Point(reach, reach, reach);
    std::vector<std::size_t> candidates;
    const auto lo = cell_of(box.min()), hi = cell_of(box.max());
    for (long i = lo[0]; i <= hi[0]; ++i)
      for (long j = lo[1]; j <= hi[1]; ++j)
        for (long k = lo[2]; k <= hi[2]; ++k)
          if (const auto it = bins.find({i, j, k}); it != bins.end())
            for (const auto f : it->second)
              if (box.intersects(faces[f].box))
                candidates.push_back(f);
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    faces_host.insert(faces_host.end(), candidates.begin(), candidates.end());
    offsets_host[cid + 1] = faces_host.size();
  }

  // To device
  vertices = ::Kokkos::View<Real * [3][3], ::Kokkos::LayoutRight>("dem_wall_vertices", n);
  num_vertices = ::Kokkos::View<unsigned int *>("dem_wall_num_vertices", n);
  normals = ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight>("dem_wall_normals", n);
  boundary = ::Kokkos::View<std::size_t *>("dem_wall_boundary", n);
  elem_offsets = ::Kokkos::View<std::size_t *>("dem_wall_elem_offsets", num_elems + 1);
  elem_faces = ::Kokkos::View<std::size_t *>("dem_wall_elem_faces", faces_host.size());
  auto vertices_h = ::Kokkos::create_mirror_view(vertices);
  auto num_vertices_h = ::Kokkos::create_mirror_view(num_vertices);
  auto normals_h = ::Kokkos::create_mirror_view(normals);
  auto boundary_h = ::Kokkos::create_mirror_view(boundary);
  auto elem_offsets_h = ::Kokkos::create_mirror_view(elem_offsets);
  auto elem_faces_h = ::Kokkos::create_mirror_view(elem_faces);
  for (const auto f : index_range(faces))
  {
    for (const auto k : make_range(3))
      for (const auto c : make_range(3))
        vertices_h(f, k, c) = faces[f].v[k](c);
    num_vertices_h(f) = faces[f].nv;
    for (const auto c : make_range(3))
      normals_h(f, c) = faces[f].normal(c);
    boundary_h(f) = faces[f].boundary;
  }
  for (const auto e : index_range(offsets_host))
    elem_offsets_h(e) = offsets_host[e];
  for (const auto k : index_range(faces_host))
    elem_faces_h(k) = faces_host[k];
  ::Kokkos::deep_copy(vertices, vertices_h);
  ::Kokkos::deep_copy(num_vertices, num_vertices_h);
  ::Kokkos::deep_copy(normals, normals_h);
  ::Kokkos::deep_copy(boundary, boundary_h);
  ::Kokkos::deep_copy(elem_offsets, elem_offsets_h);
  ::Kokkos::deep_copy(elem_faces, elem_faces_h);
}

} // namespace DEM
