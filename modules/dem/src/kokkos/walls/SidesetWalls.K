#include "SidesetWalls.h"

#include "MooseMesh.h"

#include "libmesh/elem.h"
#include "timpi/parallel_sync.h"

#include <map>
#include <algorithm>

namespace DEM
{

void
SidesetWalls::build(const MooseMesh & mesh,
                    const std::vector<const libMesh::Elem *> & cid_to_elem,
                    const std::size_t num_local,
                    const std::vector<BoundaryID> & boundaries,
                    const Real reach,
                    const Real curvature,
                    const std::vector<libMesh::BoundingBox> & rank_boxes,
                    const libMesh::Parallel::Communicator & comm)
{
  const auto num_elems = cid_to_elem.size();

  // Faces on host: vertices, vertex count, inward normal, global ID, nodes, bounding box; and as
  // a flat record of Reals for the exchange (the IDs are exact below 2^53)
  struct Face
  {
    Point v[3];
    unsigned int nv;
    Point normal;
    int64_t id;
    std::array<libMesh::dof_id_type, 3> nodes;
    libMesh::BoundingBox box;
  };
  constexpr std::size_t record_size = 9 + 1 + 3 + 1 + 3;
  const auto pack = [](const Face & face, std::vector<Real> & out)
  {
    for (const auto k : make_range(3))
      for (const auto c : make_range(3))
        out.push_back(face.v[k](c));
    out.push_back(face.nv);
    for (const auto c : make_range(3))
      out.push_back(face.normal(c));
    out.push_back(face.id);
    for (const auto k : make_range(3))
      out.push_back(face.nodes[k]);
  };
  const auto unpack = [](const Real * record)
  {
    Face face;
    for (const auto k : make_range(3))
      for (const auto c : make_range(3))
        face.v[k](c) = record[3 * k + c];
    face.nv = record[9];
    for (const auto c : make_range(3))
      face.normal(c) = record[10 + c];
    face.id = record[13];
    for (const auto k : make_range(3))
      face.nodes[k] = record[14 + k];
    for (const auto k : make_range(face.nv))
      face.box.union_with(face.v[k]);
    return face;
  };

  // The faces of this rank's own elements
  std::vector<Face> mine;
  const auto add = [&mine](const Point & a,
                           const Point & b,
                           const Point & c,
                           const unsigned int nv,
                           const Point & centroid,
                           const int64_t id,
                           const std::array<libMesh::dof_id_type, 3> & nodes)
  {
    Face face{{a, b, c}, nv, Point(), id, nodes, libMesh::BoundingBox()};
    // Normal toward the element's centroid, which is inside the domain
    face.normal = nv == 3 ? (b - a).cross(c - a) : Point(-(b - a)(1), (b - a)(0), 0);
    face.normal /= face.normal.norm();
    if (face.normal * (centroid - a) < 0)
      face.normal *= -1;
    for (const auto k : make_range(nv))
      face.box.union_with(face.v[k]);
    mine.push_back(face);
  };
  for (const auto cid : make_range(num_local))
  {
    const auto * elem = cid_to_elem[cid];
    for (const auto s : make_range(elem->n_sides()))
    {
      const auto ids = mesh.getBoundaryIDs(elem, s);
      if (std::none_of(ids.begin(), ids.end(), [&boundaries](auto id) {
            return std::find(boundaries.begin(), boundaries.end(), id) != boundaries.end();
          }))
        continue;
      const auto side = elem->side_ptr(s);
      const auto centroid = elem->vertex_average();
      const int64_t base = (static_cast<int64_t>(elem->id()) * 8 + s) * 2;
      // Vertices only: higher-order sides are used through their corner nodes
      std::vector<Point> v;
      std::vector<libMesh::dof_id_type> nd;
      for (const auto k : make_range(side->n_vertices()))
      {
        v.push_back(side->point(k));
        nd.push_back(side->node_id(k));
      }
      if (v.size() == 2)
        add(v[0], v[1], v[1], 2, centroid, base, {nd[0], nd[1], nd[1]});
      else if (v.size() == 3)
        add(v[0], v[1], v[2], 3, centroid, base, {nd[0], nd[1], nd[2]});
      else if (v.size() == 4)
      {
        add(v[0], v[1], v[2], 3, centroid, base, {nd[0], nd[1], nd[2]});
        add(v[0], v[2], v[3], 3, centroid, base + 1, {nd[0], nd[2], nd[3]});
      }
      else
        mooseError("Unsupported side with ", v.size(), " vertices for a sideset wall");
    }
  }

  // Every rank gets the faces within the reach of its bounding box, this rank's own included:
  // a collective exchange every rank enters
  std::map<processor_id_type, std::vector<Real>> to_send;
  for (const auto & face : mine)
  {
    auto box = face.box;
    box.first -= Point(reach, reach, reach);
    box.second += Point(reach, reach, reach);
    for (const auto rank : index_range(rank_boxes))
      if (box.intersects(rank_boxes[rank]))
        pack(face, to_send[rank]);
  }
  std::vector<Face> faces;
  Parallel::push_parallel_vector_data(
      comm,
      to_send,
      [&faces, &unpack](processor_id_type, const std::vector<Real> & data)
      {
        for (std::size_t r = 0; r < data.size() / record_size; ++r)
          faces.push_back(unpack(data.data() + r * record_size));
      });
  // Deterministic order whatever the arrival order
  std::sort(faces.begin(), faces.end(), [](const Face & a, const Face & b) { return a.id < b.id; });
  n = faces.size();
  nodes.resize(n);
  for (const auto f : index_range(faces))
    nodes[f] = faces[f].nodes;

  // Adjacency across edges, by the node pair of each edge (the faces of a split quad share its
  // diagonal); segments in 2D are adjacent through their end nodes. Two neighbors are one
  // surface when their normals are within the curvature angle
  std::vector<std::array<int, 3>> neighbors_host(n, {-1, -1, -1});
  std::vector<std::array<bool, 3>> same_host(n, {false, false, false});
  std::map<std::pair<libMesh::dof_id_type, libMesh::dof_id_type>, std::vector<std::pair<std::size_t, unsigned int>>> edges;
  for (const auto f : index_range(faces))
    if (faces[f].nv == 3)
      for (const auto e : make_range(3u))
      {
        const auto a = faces[f].nodes[e], b = faces[f].nodes[(e + 1) % 3];
        edges[{std::min(a, b), std::max(a, b)}].emplace_back(f, e);
      }
    else
      for (const auto e : make_range(2u))
        edges[{faces[f].nodes[e], faces[f].nodes[e]}].emplace_back(f, e);
  const Real cos_curvature = std::cos(curvature);
  for (const auto & [edge, sharing] : edges)
    if (sharing.size() == 2)
    {
      const auto [f, e] = sharing[0];
      const auto [g, d] = sharing[1];
      neighbors_host[f][e] = g;
      neighbors_host[g][d] = f;
      const bool same = faces[f].normal * faces[g].normal >= cos_curvature;
      same_host[f][e] = same;
      same_host[g][d] = same;
    }

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
  velocities = ::Kokkos::View<Real * [3][3], ::Kokkos::LayoutRight>("dem_wall_velocities", n);
  overflow = ::Kokkos::View<int>("dem_wall_overflow");
  moving = false;
  ids = ::Kokkos::View<int64_t *>("dem_wall_ids", n);
  neighbors = ::Kokkos::View<int * [3], ::Kokkos::LayoutRight>("dem_wall_neighbors", n);
  same_surface = ::Kokkos::View<bool * [3], ::Kokkos::LayoutRight>("dem_wall_same_surface", n);
  elem_offsets = ::Kokkos::View<std::size_t *>("dem_wall_elem_offsets", num_elems + 1);
  elem_faces = ::Kokkos::View<std::size_t *>("dem_wall_elem_faces", faces_host.size());
  auto vertices_h = ::Kokkos::create_mirror_view(vertices);
  auto num_vertices_h = ::Kokkos::create_mirror_view(num_vertices);
  auto normals_h = ::Kokkos::create_mirror_view(normals);
  auto ids_h = ::Kokkos::create_mirror_view(ids);
  auto neighbors_h = ::Kokkos::create_mirror_view(neighbors);
  auto same_h = ::Kokkos::create_mirror_view(same_surface);
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
    ids_h(f) = faces[f].id;
    for (const auto e : make_range(3))
    {
      neighbors_h(f, e) = neighbors_host[f][e];
      same_h(f, e) = same_host[f][e];
    }
  }
  for (const auto e : index_range(offsets_host))
    elem_offsets_h(e) = offsets_host[e];
  for (const auto k : index_range(faces_host))
    elem_faces_h(k) = faces_host[k];
  ::Kokkos::deep_copy(vertices, vertices_h);
  ::Kokkos::deep_copy(num_vertices, num_vertices_h);
  ::Kokkos::deep_copy(normals, normals_h);
  ::Kokkos::deep_copy(ids, ids_h);
  ::Kokkos::deep_copy(neighbors, neighbors_h);
  ::Kokkos::deep_copy(same_surface, same_h);
  ::Kokkos::deep_copy(elem_offsets, elem_offsets_h);
  ::Kokkos::deep_copy(elem_faces, elem_faces_h);
}

void
SidesetWalls::move(const std::vector<Point> & positions, const Real dt)
{
  mooseAssert(positions.size() == 3 * n, "One position per face vertex");
  auto vertices_h = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, vertices);
  auto velocities_h = ::Kokkos::create_mirror_view(velocities);
  auto normals_h = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, normals);
  auto num_vertices_h = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, num_vertices);
  moving = false;
  for (const auto f : make_range(n))
  {
    for (const auto k : make_range(3))
      for (const auto c : make_range(3))
      {
        velocities_h(f, k, c) = dt > 0 ? (positions[3 * f + k](c) - vertices_h(f, k, c)) / dt : 0;
        if (dt <= 0)
          vertices_h(f, k, c) = positions[3 * f + k](c);
        moving = moving || velocities_h(f, k, c) != 0;
      }
    // The normal of the face at its new position, oriented like the current one
    const Point & a = positions[3 * f], & b = positions[3 * f + 1], & c = positions[3 * f + 2];
    Point normal = num_vertices_h(f) == 2 ? Point(-(b - a)(1), (b - a)(0), 0) : (b - a).cross(c - a);
    if (normal.norm() > 0)
    {
      normal /= normal.norm();
      if (normal * Point(normals_h(f, 0), normals_h(f, 1), normals_h(f, 2)) < 0)
        normal *= -1;
      for (const auto c : make_range(3))
        normals_h(f, c) = normal(c);
    }
  }
  ::Kokkos::deep_copy(velocities, velocities_h);
  ::Kokkos::deep_copy(normals, normals_h);
  if (dt <= 0)
    ::Kokkos::deep_copy(vertices, vertices_h);
}

void
SidesetWalls::save(std::vector<Real> & positions, std::vector<Real> & velocities_out) const
{
  auto vertices_h = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, vertices);
  auto velocities_h = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, velocities);
  positions.assign(vertices_h.data(), vertices_h.data() + 9 * n);
  velocities_out.assign(velocities_h.data(), velocities_h.data() + 9 * n);
}

void
SidesetWalls::load(const std::vector<Real> & positions, const std::vector<Real> & velocities_in)
{
  mooseAssert(positions.size() == 9 * n && velocities_in.size() == 9 * n,
              "Checkpoint of a different wall mesh");
  auto vertices_h = ::Kokkos::create_mirror_view(vertices);
  auto velocities_h = ::Kokkos::create_mirror_view(velocities);
  moving = false;
  for (const auto f : make_range(n))
    for (const auto k : make_range(3))
      for (const auto c : make_range(3))
      {
        vertices_h(f, k, c) = positions[9 * f + 3 * k + c];
        velocities_h(f, k, c) = velocities_in[9 * f + 3 * k + c];
        moving = moving || velocities_h(f, k, c) != 0;
      }
  ::Kokkos::deep_copy(vertices, vertices_h);
  ::Kokkos::deep_copy(velocities, velocities_h);
}

void
SidesetWalls::advance(const Real dt)
{
  if (!moving)
    return;
  const auto v = vertices;
  const auto vel = velocities;
  ::Kokkos::parallel_for(
      "dem_wall_advance", n, KOKKOS_LAMBDA(const std::size_t f) {
        for (unsigned int k = 0; k < 3; ++k)
          for (unsigned int c = 0; c < 3; ++c)
            v(f, k, c) += dt * vel(f, k, c);
      });
}

} // namespace DEM
