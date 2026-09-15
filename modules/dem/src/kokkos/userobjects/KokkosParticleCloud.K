#include "KokkosParticleCloud.h"

#include "MooseMesh.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h"
#include "timpi/parallel_sync.h"

#include <algorithm>

registerMooseObject("DemApp", KokkosParticleCloud);

InputParameters
KokkosParticleCloud::validParams()
{
  InputParameters params = Moose::Kokkos::GeneralUserObject::validParams();
  params.addClassDescription(
      "Integrates a cloud of spherical particles on device with Velocity Verlet over a number of "
      "substeps per time step, tracking the element containing each particle by walking across "
      "element faces and migrating particles between ranks as they cross partitions.");
  params.addRequiredParam<std::vector<Point>>("initial_positions",
                                              "Initial positions of all particles");
  params.addParam<RealVectorValue>(
      "initial_velocity", RealVectorValue(0, 0, 0), "Initial velocity of every particle");
  params.addParam<std::vector<Point>>(
      "initial_velocities",
      "Per-particle initial velocities overriding initial_velocity, one per initial position");
  params.addParam<RealVectorValue>("initial_angular_velocity",
                                   RealVectorValue(0, 0, 0),
                                   "Initial body-frame angular velocity of every particle");
  params.addRequiredRangeCheckedParam<Real>("radius", "radius > 0", "Particle radius");
  params.addParam<std::vector<Real>>(
      "initial_radii", "Per-particle radii overriding radius, one per initial position");
  params.addRequiredRangeCheckedParam<Real>("density", "density > 0", "Particle density");
  params.addParam<RealVectorValue>(
      "gravity", RealVectorValue(0, 0, 0), "Gravitational acceleration");
  params.addRangeCheckedParam<Real>("normal_stiffness",
                                    0,
                                    "normal_stiffness >= 0",
                                    "Spring stiffness of the linear spring-dashpot normal contact; "
                                    "zero disables contact");
  params.addRangeCheckedParam<Real>("normal_damping",
                                    0,
                                    "normal_damping >= 0",
                                    "Dashpot coefficient of the linear spring-dashpot normal contact");
  params.addParam<std::vector<Point>>(
      "wall_points", {}, "A point on each fixed planar wall, one per entry of wall_normals");
  params.addParam<std::vector<RealVectorValue>>(
      "wall_normals",
      {},
      "Normal of each fixed planar wall pointing into the domain, one per entry of wall_points; "
      "particles overlapping a wall get the normal contact force against it");
  params.addRangeCheckedParam<unsigned int>(
      "substeps", 1, "substeps > 0", "Number of integration substeps per time step");
  // A particle should cross only a few cells per substep; a much larger walk indicates a mesh the
  // face-plane test cannot handle, so bound it rather than loop for the mesh size
  params.addParam<unsigned int>(
      "max_hops", 32, "Maximum number of face-neighbor hops per particle per substep");
  params.addRangeCheckedParam<Real>("skin",
                                    "skin >= 0",
                                    "Verlet skin distance of the neighbor list; the list is rebuilt "
                                    "once some particle has moved skin / 2. Defaults to a tenth of "
                                    "radius");
  params.addParam<bool>("gpu_aware_mpi",
                        false,
                        "Pass device buffers directly to MPI for the ghost update on a GPU build; "
                        "requires a GPU-aware MPI library");
  params.addParam<bool>("verify",
                        false,
                        "Check the device element assignment against libMesh's PointLocator "
                        "every step and error on any mismatch");
  params.addParam<bool>("verify_neighbor_list",
                        false,
                        "Check the neighbor list against a brute-force pair search after every "
                        "build and error on any mismatch");
  return params;
}

KokkosParticleCloud::KokkosParticleCloud(const InputParameters & parameters)
  : Moose::Kokkos::GeneralUserObject(parameters),
    _initial_positions(getParam<std::vector<Point>>("initial_positions")),
    _initial_velocity(getParam<RealVectorValue>("initial_velocity")),
    _initial_velocities(isParamValid("initial_velocities")
                            ? getParam<std::vector<Point>>("initial_velocities")
                            : std::vector<Point>()),
    _initial_angular_velocity(getParam<RealVectorValue>("initial_angular_velocity")),
    _radius(getParam<Real>("radius")),
    _initial_radii(isParamValid("initial_radii") ? getParam<std::vector<Real>>("initial_radii")
                                                 : std::vector<Real>()),
    _density(getParam<Real>("density")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _contact({getParam<Real>("normal_stiffness"), getParam<Real>("normal_damping")}),
    _walls(getParam<std::vector<Point>>("wall_points"),
           getParam<std::vector<RealVectorValue>>("wall_normals")),
    _substeps(getParam<unsigned int>("substeps")),
    _max_hops(getParam<unsigned int>("max_hops")),
    _skin(isParamValid("skin") ? getParam<Real>("skin") : 0.1 * _radius),
    _verify(getParam<bool>("verify")),
    _verify_neighbor_list(getParam<bool>("verify_neighbor_list")),
    _ghost_tag(_communicator.get_unique_tag()),
    _gpu_aware_mpi(getParam<bool>("gpu_aware_mpi"))
{
  if (_initial_radii.size() && _initial_radii.size() != _initial_positions.size())
    paramError("initial_radii", "Must have one entry per initial position");
  for (const auto r : _initial_radii)
    if (r <= 0)
      paramError("initial_radii", "Radii must be positive");
  if (_initial_velocities.size() && _initial_velocities.size() != _initial_positions.size())
    paramError("initial_velocities", "Must have one entry per initial position");
  if (getParam<std::vector<Point>>("wall_points").size() !=
      getParam<std::vector<RealVectorValue>>("wall_normals").size())
    paramError("wall_normals", "Must have one entry per wall point");
  for (const auto & normal : getParam<std::vector<RealVectorValue>>("wall_normals"))
    if (normal.norm() == 0)
      paramError("wall_normals", "Normals must be nonzero");
  if (_walls.n && _contact.stiffness == 0)
    paramError("wall_points", "Walls need a positive normal_stiffness");

  // Side centroids and normals are only cached on request
  _fe_problem.mesh().getKokkosMesh()->setNeedsElementSideGeometry();
}

void
KokkosParticleCloud::initialSetup()
{
  const auto & mesh = kokkosMesh();
  const auto num_total = mesh.getNumLocalAndPossiblyOneNeighborLayerGhostElements();

  // Reverse element ID map, owner ranks, and libMesh IDs for local and one-layer ghost elements
  _cid_to_elem.assign(num_total, nullptr);
  auto owner_host = ::Kokkos::View<processor_id_type *, ::Kokkos::HostSpace>("dem_owner", num_total);
  auto id_host = ::Kokkos::View<dof_id_type *, ::Kokkos::HostSpace>("dem_libmesh_id", num_total);
  for (const auto elem : _fe_problem.mesh().getMesh().active_local_element_ptr_range())
  {
    const auto cid = mesh.getContiguousElementID(elem);
    _cid_to_elem[cid] = elem;
    owner_host(cid) = processor_id();
    id_host(cid) = elem->id();
  }
  for (const auto & [elem, cid] : mesh.getGhostElemIdMapping())
  {
    _cid_to_elem[cid] = elem;
    owner_host(cid) = elem->processor_id();
    id_host(cid) = elem->id();
  }
  _owner = ::Kokkos::View<processor_id_type *>("dem_owner", num_total);
  _libmesh_id = ::Kokkos::View<dof_id_type *>("dem_libmesh_id", num_total);
  ::Kokkos::deep_copy(_owner, owner_host);
  ::Kokkos::deep_copy(_libmesh_id, id_host);

  // Each rank keeps the particles that start in one of its local elements
  _locator = _fe_problem.mesh().getPointLocator();
  _locator->enable_out_of_mesh_mode();

  std::vector<std::pair<int64_t, ContiguousElementID>> local;
  for (const auto i : index_range(_initial_positions))
  {
    const auto * elem = (*_locator)(_initial_positions[i]);
    if (elem && elem->processor_id() == processor_id())
      local.emplace_back(i, mesh.getContiguousElementID(elem));
  }

  const auto n = local.size();
  DEM::ParticleCloud::HostMirror host(n);
  for (const auto j : make_range(n))
  {
    const auto [gid, cid] = local[j];
    const Real radius = _initial_radii.size() ? _initial_radii[gid] : _radius;
    const Real mass = _density * 4.0 / 3.0 * libMesh::pi * Utility::pow<3>(radius);
    // Solid sphere
    const Real inertia = 0.4 * mass * radius * radius;
    host.gid(j) = gid;
    for (const auto c : make_range(3))
    {
      host.x(j, c) = _initial_positions[gid](c);
      host.v(j, c) = _initial_velocities.size() ? _initial_velocities[gid](c) : _initial_velocity(c);
      host.omega(j, c) = _initial_angular_velocity(c);
      host.f(j, c) = 0;
      host.tau(j, c) = 0;
    }
    // Identity orientation
    host.q(j, 0) = 1;
    host.q(j, 1) = host.q(j, 2) = host.q(j, 3) = 0;
    host.r(j) = radius;
    host.m(j) = mass;
    host.inertia(j) = inertia;
    host.elem(j) = cid;
    host.target_rank(j) = libMesh::DofObject::invalid_processor_id;
    host.hops(j) = 0;
  }
  _cloud.allocate(n);
  _cloud.copyFrom(host, 0);
  _cloud.n = n;
  allocateScratch(n);

  // Ghosting is decided against the neighboring ranks' bounding boxes, gathered once since the
  // mesh does not move
  _r_max = _initial_radii.size() ? *std::max_element(_initial_radii.begin(), _initial_radii.end())
                                 : _radius;
  const auto my_box = _fe_problem.mesh().getInflatedProcessorBoundingBox();
  std::vector<Real> bounds = {my_box.min()(0),
                              my_box.min()(1),
                              my_box.min()(2),
                              my_box.max()(0),
                              my_box.max()(1),
                              my_box.max()(2)};
  _communicator.allgather(bounds);
  _rank_boxes.clear();
  for (const auto rank : make_range(n_processors()))
    _rank_boxes.emplace_back(Point(bounds[6 * rank], bounds[6 * rank + 1], bounds[6 * rank + 2]),
                             Point(bounds[6 * rank + 3], bounds[6 * rank + 4], bounds[6 * rank + 5]));
  const Real cutoff = 2 * _r_max + _skin;
  auto reach = my_box;
  reach.first -= Point(cutoff, cutoff, cutoff);
  reach.second += Point(cutoff, cutoff, cutoff);
  _neighbor_ranks.clear();
  for (const auto rank : make_range(n_processors()))
    if (rank != processor_id() && reach.intersects(_rank_boxes[rank]))
      _neighbor_ranks.push_back(rank);

  updateNeighborList();
  computeForces();
  // Reduce so the initial values are global before the first execution
  count();
  finalize();
}

void
KokkosParticleCloud::compute()
{
  const auto cloud = _cloud;
  ::Kokkos::parallel_for(
      "dem_reset_hops", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) { cloud.hops(i) = 0; });

  // Velocity Verlet; the forces at the start of the step were computed at the end of the last one
  for (unsigned int substep = 0; substep < _substeps; ++substep)
  {
    kickDrift();
    walk();
    if (!updateNeighborList())
      startGhostUpdate();
    computeForces();
    kick();
  }

  // Particles the walk could not resolve get a second chance from the point locator, and may
  // turn out to belong to another rank; the walk that follows flags those for migration
  resolveUnresolved();
  walk();

  // Remove the exited particles and hand the ones that reached a ghost element to their owner,
  // who continues the walk (with the same fallback, since the hop budget may not reach from the
  // element the sender saw), until every particle has settled. Each round moves a particle at
  // least one partition, so the number of rounds is bounded by the number of ranks
  _num_migrated = 0;
  for (const auto round : make_range(n_processors()))
  {
    auto n_sent = migrate();
    _num_migrated += n_sent;
    _communicator.sum(n_sent);
    if (!n_sent)
      break;
    walk();
    resolveUnresolved();
    walk();

    if (round + 1 == n_processors())
      mooseError("Particles were still crossing partitions after ", n_processors(), " rounds");
  }
  // Received particles arrive without forces, and the next step's first half-kick needs them
  if (!updateNeighborList())
    startGhostUpdate();
  computeForces();

  count();
  if (_verify)
    verify();
}

void
KokkosParticleCloud::computeForces(const ::Kokkos::View<std::size_t *> & subset,
                                   const std::size_t count)
{
  const auto cloud = _cloud;
  const Moose::Kokkos::Real3 g(_gravity(0), _gravity(1), _gravity(2));
  const auto contact = _contact;
  const auto walls = _walls;
  const auto offsets = _neighbor_list.offsets();
  const auto pairs = _neighbor_list.pairs();
  constexpr auto exited = DEM::ParticleCloud::exited;

  ::Kokkos::parallel_for(
      "dem_forces", count, KOKKOS_LAMBDA(const std::size_t s) {
        const auto i = subset(s);
        Moose::Kokkos::Real3 f = 0;
        // Exited particles are inert until they are removed
        if (cloud.elem(i) != exited)
        {
          f = cloud.m(i) * g;

          // Sphere-sphere narrow phase over the neighbor list; every pair is listed from both
          // sides so each particle sums its own contacts
          if (contact.stiffness > 0)
          {
            const Moose::Kokkos::Real3 xi(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2));
            const Moose::Kokkos::Real3 vi(cloud.v(i, 0), cloud.v(i, 1), cloud.v(i, 2));
            for (auto k = offsets(i); k < offsets(i + 1); ++k)
            {
              const auto j = pairs(k);
              if (cloud.elem(j) == exited)
                continue;
              const Moose::Kokkos::Real3 xj(cloud.x(j, 0), cloud.x(j, 1), cloud.x(j, 2));
              const Moose::Kokkos::Real3 d = xi - xj;
              const Real distance = d.norm();
              const Real overlap = cloud.r(i) + cloud.r(j) - distance;
              if (overlap <= 0)
                continue;
              const Moose::Kokkos::Real3 normal = (1.0 / distance) * d;
              const Moose::Kokkos::Real3 vj(cloud.v(j, 0), cloud.v(j, 1), cloud.v(j, 2));
              const Real normal_velocity = (vi - vj).dot_product(normal);
              f += contact.normalForce(overlap, normal_velocity) * normal;
            }

            // Sphere-plane narrow phase against every wall, which is at rest
            for (std::size_t w = 0; w < walls.n; ++w)
            {
              const Real overlap = walls.overlap(w, xi, cloud.r(i));
              if (overlap <= 0)
                continue;
              const auto normal = walls.normal(w);
              f += contact.normalForce(overlap, vi.dot_product(normal)) * normal;
            }
          }
        }

        for (unsigned int c = 0; c < 3; ++c)
        {
          cloud.f(i, c) = f(c);
          cloud.tau(i, c) = 0;
        }
      });
}

void
KokkosParticleCloud::computeForces()
{
  // Interior particles do not need the ghost update, so their forces overlap its communication
  computeForces(_neighbor_list.interior(), _neighbor_list.numInterior());
  finishGhostUpdate();
  computeForces(_neighbor_list.boundary(), _neighbor_list.numBoundary());
}

void
KokkosParticleCloud::kickDrift()
{
  const auto cloud = _cloud;
  const Real dt = _dt / _substeps;
  constexpr auto exited = DEM::ParticleCloud::exited;

  ::Kokkos::parallel_for(
      "dem_kick_drift", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.elem(i) == exited)
          return;
        const Real half_dt_over_m = 0.5 * dt / cloud.m(i);
        const Real half_dt_over_I = 0.5 * dt / cloud.inertia(i);
        auto x = libMesh::Kokkos::make_vector_ref(cloud.x, i);
        auto v = libMesh::Kokkos::make_vector_ref(cloud.v, i);
        auto omega = libMesh::Kokkos::make_vector_ref(cloud.omega, i);

        v.add_scaled(libMesh::Kokkos::make_vector_ref(cloud.f, i), half_dt_over_m);
        x.add_scaled(v, dt);
        omega.add_scaled(libMesh::Kokkos::make_vector_ref(cloud.tau, i), half_dt_over_I);

        // Rotate the orientation by the body-frame angular velocity over the substep, exactly for
        // constant angular velocity: q <- q * (cos(theta/2), sin(theta/2) * omega_hat)
        const Real omega_norm = omega.norm();
        if (omega_norm > 0)
        {
          const Real half_theta = 0.5 * omega_norm * dt;
          const Real s = ::Kokkos::sin(half_theta) / omega_norm;
          const Real dw = ::Kokkos::cos(half_theta);
          const Real dx = s * omega(0), dy = s * omega(1), dz = s * omega(2);
          const Real qw = cloud.q(i, 0), qx = cloud.q(i, 1), qy = cloud.q(i, 2),
                     qz = cloud.q(i, 3);
          cloud.q(i, 0) = qw * dw - qx * dx - qy * dy - qz * dz;
          cloud.q(i, 1) = qw * dx + qx * dw + qy * dz - qz * dy;
          cloud.q(i, 2) = qw * dy - qx * dz + qy * dw + qz * dx;
          cloud.q(i, 3) = qw * dz + qx * dy - qy * dx + qz * dw;
          // Renormalize to keep roundoff from accumulating in the unit quaternion
          const Real inv_norm =
              1.0 / ::Kokkos::sqrt(cloud.q(i, 0) * cloud.q(i, 0) + cloud.q(i, 1) * cloud.q(i, 1) +
                                   cloud.q(i, 2) * cloud.q(i, 2) + cloud.q(i, 3) * cloud.q(i, 3));
          for (unsigned int c = 0; c < 4; ++c)
            cloud.q(i, c) *= inv_norm;
        }
      });
}

void
KokkosParticleCloud::kick()
{
  const auto cloud = _cloud;
  const Real dt = _dt / _substeps;
  constexpr auto exited = DEM::ParticleCloud::exited;

  ::Kokkos::parallel_for(
      "dem_kick", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.elem(i) == exited)
          return;
        auto v = libMesh::Kokkos::make_vector_ref(cloud.v, i);
        auto omega = libMesh::Kokkos::make_vector_ref(cloud.omega, i);
        v.add_scaled(libMesh::Kokkos::make_vector_ref(cloud.f, i), 0.5 * dt / cloud.m(i));
        omega.add_scaled(libMesh::Kokkos::make_vector_ref(cloud.tau, i),
                         0.5 * dt / cloud.inertia(i));
      });
}

void
KokkosParticleCloud::walk()
{
  // Shallow copies captured by value into the kernel
  const Moose::Kokkos::Mesh mesh = kokkosMesh();
  const auto num_local = mesh.getNumLocalElements();
  const auto cloud = _cloud;
  const auto owner = _owner;
  const auto max_hops = _max_hops;
  constexpr auto invalid_elem = libMesh::DofObject::invalid_id;
  constexpr auto unresolved = DEM::ParticleCloud::unresolved;
  constexpr auto exited = DEM::ParticleCloud::exited;
  constexpr auto invalid_proc = libMesh::DofObject::invalid_processor_id;

  ::Kokkos::parallel_for(
      "dem_walk", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        auto elem = cloud.elem(i);
        cloud.target_rank(i) = invalid_proc;
        if (elem == unresolved || elem == exited)
          return;

        const Moose::Kokkos::Real3 xp(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2));
        unsigned int hops = 0;
        while (true)
        {
          // Ghost elements have no cached side geometry; the owner continues the walk after
          // migration
          if (elem >= num_local)
          {
            cloud.target_rank(i) = owner(elem);
            break;
          }

          // The particle is inside if it is behind every outward face plane; otherwise leave
          // through the face it is furthest in front of
          const auto & info = mesh.getElementInfo(elem);
          Real d_max = 0;
          unsigned int exit_side = libMesh::invalid_uint;
          // Raw loop: make_range() is not device-callable on GPU backends
          for (unsigned int side = 0; side < mesh.getNumSides(info.type); ++side)
          {
            const Real d =
                (xp - mesh.getSideCentroid(elem, side)).dot_product(mesh.getSideNormal(elem, side));
            if (d > d_max)
            {
              d_max = d;
              exit_side = side;
            }
          }
          if (exit_side == libMesh::invalid_uint)
            break;

          if (hops == max_hops)
          {
            elem = unresolved;
            break;
          }
          elem = mesh.getNeighbor(elem, exit_side);
          ++hops;
          if (elem == invalid_elem)
          {
            // Left through a face without a neighbor, i.e. the mesh boundary
            elem = exited;
            break;
          }
        }
        cloud.elem(i) = elem;
        cloud.hops(i) += hops;
      });
}

void
KokkosParticleCloud::resolveUnresolved()
{
  const auto cloud = _cloud;
  constexpr auto unresolved = DEM::ParticleCloud::unresolved;
  std::size_t n_unresolved = 0;
  ::Kokkos::parallel_reduce(
      "dem_count_unresolved",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum) {
        sum += cloud.elem(i) == unresolved;
      },
      n_unresolved);
  if (!n_unresolved)
    return;

  // Only the unresolved particles cross to the host; a particle found in a ghost element is
  // handed off by the walk that follows, one found beyond the ghost layer stays unresolved
  const auto & mesh = kokkosMesh();
  const auto & ghost_ids = mesh.getGhostElemIdMapping();
  auto x_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.x);
  auto elem_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.elem);
  for (const auto i : make_range(_cloud.size()))
  {
    if (elem_host(i) != unresolved)
      continue;
    const auto * elem = (*_locator)(Point(x_host(i, 0), x_host(i, 1), x_host(i, 2)));
    if (!elem)
      elem_host(i) = DEM::ParticleCloud::exited;
    else if (elem->processor_id() == processor_id())
      elem_host(i) = mesh.getContiguousElementID(elem);
    else if (const auto it = ghost_ids.find(elem); it != ghost_ids.end())
      elem_host(i) = it->second;
  }
  ::Kokkos::deep_copy(_cloud.elem, elem_host);
}

bool
KokkosParticleCloud::updateNeighborList()
{
  // Both sides of a cross-rank pair must rebuild together
  unsigned int stale = _neighbor_list.stale(_cloud, _cloud.size());
  if (n_processors() > 1)
    _communicator.max(stale);
  if (!stale)
    return false;

  // A rebuild replaces the ghosts, so a pending update must land first
  finishGhostUpdate();
  exchangeGhosts();

  // Ghosts lie within the cutoff of our box; particles that walked out of the mesh are clamped
  // into the edge bins
  const Real cutoff = 2 * _r_max + _skin;
  auto box = _rank_boxes[processor_id()];
  box.first -= Point(cutoff, cutoff, cutoff);
  box.second += Point(cutoff, cutoff, cutoff);
  _neighbor_list.build(_cloud, _cloud.size(), _cloud.size() + _num_ghosts, box.min(), box.max(), _skin);
  ++_num_neighbor_list_builds;
  if (_verify_neighbor_list)
    verifyNeighborList();
  return true;
}

void
KokkosParticleCloud::exchangeGhosts()
{
  constexpr auto record_size = DEM::ParticleCloud::record_size;
  constexpr auto exited = DEM::ParticleCloud::exited;
  const auto cloud = _cloud;
  const auto n = _cloud.size();
  const Real cutoff = 2 * _r_max + _skin;
  _num_ghosts = 0;
  _ghost_send_counts.assign(_neighbor_ranks.size(), 0);
  _ghost_recv.clear();
  if (_neighbor_ranks.empty())
    return;

  // Flag, per neighbor rank, the local particles within the cutoff of that rank's box, then
  // compact the flags into the send index with a scan
  const auto num_neighbors = _neighbor_ranks.size();
  if (_slot.extent(0) < n * num_neighbors)
    _slot = ::Kokkos::View<std::size_t *>("dem_migrate_slot", n * num_neighbors);
  const auto flags = _slot;
  ::Kokkos::View<Real * [6], ::Kokkos::LayoutRight> boxes("dem_ghost_boxes", num_neighbors);
  auto boxes_host = ::Kokkos::create_mirror_view(boxes);
  for (const auto b : index_range(_neighbor_ranks))
    for (const auto c : make_range(3))
    {
      boxes_host(b, c) = _rank_boxes[_neighbor_ranks[b]].min()(c) - cutoff;
      boxes_host(b, 3 + c) = _rank_boxes[_neighbor_ranks[b]].max()(c) + cutoff;
    }
  ::Kokkos::deep_copy(boxes, boxes_host);
  ::Kokkos::parallel_for(
      "dem_ghost_flag", n * num_neighbors, KOKKOS_LAMBDA(const std::size_t k) {
        const auto b = k / n, i = k % n;
        bool inside = cloud.elem(i) != exited;
        for (unsigned int c = 0; c < 3; ++c)
          inside = inside && cloud.x(i, c) >= boxes(b, c) && cloud.x(i, c) <= boxes(b, 3 + c);
        flags(k) = inside;
      });
  std::size_t total = 0;
  ::Kokkos::parallel_scan(
      "dem_ghost_scan",
      n * num_neighbors,
      KOKKOS_LAMBDA(const std::size_t k, std::size_t & sum, const bool final) {
        const auto flag = flags(k);
        if (final)
          flags(k) = flag ? sum : std::numeric_limits<std::size_t>::max();
        sum += flag;
      },
      total);
  if (_ghost_send_index.extent(0) < total)
    _ghost_send_index = ::Kokkos::View<std::size_t *>("dem_ghost_send_index", total);
  if (_ghost_buffer.extent(0) < total * record_size)
  {
    _ghost_buffer = ::Kokkos::View<Real *>("dem_ghost_buffer", total * record_size);
    _ghost_buffer_host = ::Kokkos::create_mirror_view(_ghost_buffer);
  }
  const auto send_index = _ghost_send_index;
  const auto buffer = _ghost_buffer;
  ::Kokkos::parallel_for(
      "dem_ghost_pack", n * num_neighbors, KOKKOS_LAMBDA(const std::size_t k) {
        const auto slot = flags(k);
        if (slot == std::numeric_limits<std::size_t>::max())
          return;
        const auto i = k % n;
        send_index(slot) = i;
        cloud.pack(i, 0, &buffer(slot * record_size));
      });

  // The send index is ordered by neighbor rank, so the per-rank counts are the slot boundaries
  auto buffer_host = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{},
      ::Kokkos::subview(buffer, std::make_pair(std::size_t(0), total * record_size)));
  auto flags_host = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{}, ::Kokkos::subview(flags, std::make_pair(std::size_t(0), n * num_neighbors)));
  std::map<processor_id_type, std::vector<Real>> to_send;
  for (const auto b : index_range(_neighbor_ranks))
  {
    auto & data = to_send[_neighbor_ranks[b]];
    for (const auto i : make_range(n))
    {
      const auto slot = flags_host(b * n + i);
      if (slot == std::numeric_limits<std::size_t>::max())
        continue;
      data.insert(data.end(), &buffer_host(slot * record_size), &buffer_host(slot * record_size) + record_size);
      ++_ghost_send_counts[b];
    }
  }

  std::map<processor_id_type, std::vector<Real>> received;
  Parallel::push_parallel_vector_data(
      _communicator,
      to_send,
      [&received](processor_id_type source, const std::vector<Real> & data)
      { received[source] = data; });

  // Append the ghosts behind the local particles, grouped by source rank in rank order so the
  // forwarding of each substep lands in fixed slots
  for (const auto & [source, data] : received)
    _num_ghosts += data.size() / record_size;
  const auto end = n + _num_ghosts;
  if (end > _cloud.capacity())
    reserve(std::max(end, 2 * _cloud.capacity()));
  std::size_t begin = n;
  for (const auto & [source, data] : received)
  {
    const auto count = data.size() / record_size;
    _ghost_recv.push_back({source, begin, count});
    DEM::ParticleCloud::HostMirror host(count);
    for (const auto r : make_range(count))
    {
      host.unpack(data.data() + r * record_size, r);
      host.elem(r) = DEM::ParticleCloud::unresolved;
    }
    _cloud.copyFrom(host, begin);
    begin += count;
  }
  if (_ghost_recv_buffer.extent(0) < _num_ghosts * 6)
  {
    _ghost_recv_buffer = ::Kokkos::View<Real *>("dem_ghost_recv_buffer", _num_ghosts * 6);
    _ghost_recv_buffer_host = ::Kokkos::create_mirror_view(_ghost_recv_buffer);
  }
}

void
KokkosParticleCloud::startGhostUpdate()
{
  if (_neighbor_ranks.empty())
    return;
  ++_num_ghost_forwards;

  // Positions and velocities only; the rest of a ghost's state is fixed between rebuilds
  constexpr std::size_t stride = 6;
  const auto cloud = _cloud;
  const auto send_index = _ghost_send_index;
  const auto buffer = _ghost_buffer;
  std::size_t total = 0;
  for (const auto count : _ghost_send_counts)
    total += count;
  ::Kokkos::parallel_for(
      "dem_ghost_forward_pack", total, KOKKOS_LAMBDA(const std::size_t slot) {
        const auto i = send_index(slot);
        for (unsigned int c = 0; c < 3; ++c)
        {
          buffer(slot * stride + c) = cloud.x(i, c);
          buffer(slot * stride + 3 + c) = cloud.v(i, c);
        }
      });
  ::Kokkos::fence();

  // MPI reads the device buffers directly when the host can access them or the MPI library is
  // GPU-aware, and host mirrors otherwise
  const bool direct = _device_buffers_to_mpi || _gpu_aware_mpi;
  Real * send_ptr = _ghost_buffer.data();
  Real * recv_ptr = _ghost_recv_buffer.data();
  if (!direct)
  {
    ::Kokkos::deep_copy(::Kokkos::subview(_ghost_buffer_host, std::make_pair(std::size_t(0), total * stride)),
                        ::Kokkos::subview(_ghost_buffer, std::make_pair(std::size_t(0), total * stride)));
    send_ptr = _ghost_buffer_host.data();
    recv_ptr = _ghost_recv_buffer_host.data();
  }

  _ghost_requests.clear();
  const auto n = _cloud.size();
  for (const auto & block : _ghost_recv)
  {
    _ghost_requests.emplace_back();
    MPI_Irecv(recv_ptr + (block.begin - n) * stride,
              block.count * stride,
              MPI_DOUBLE,
              block.source,
              _ghost_tag.value(),
              _communicator.get(),
              &_ghost_requests.back());
  }
  std::size_t slot = 0;
  for (const auto b : index_range(_neighbor_ranks))
  {
    if (_ghost_send_counts[b])
    {
      _ghost_requests.emplace_back();
      MPI_Isend(send_ptr + slot * stride,
                _ghost_send_counts[b] * stride,
                MPI_DOUBLE,
                _neighbor_ranks[b],
                _ghost_tag.value(),
                _communicator.get(),
                &_ghost_requests.back());
    }
    slot += _ghost_send_counts[b];
  }
}

void
KokkosParticleCloud::finishGhostUpdate()
{
  if (_ghost_requests.empty())
    return;
  MPI_Waitall(_ghost_requests.size(), _ghost_requests.data(), MPI_STATUSES_IGNORE);
  _ghost_requests.clear();

  constexpr std::size_t stride = 6;
  const auto n = _cloud.size();
  const auto num_ghosts = _num_ghosts;
  if (!(_device_buffers_to_mpi || _gpu_aware_mpi))
    ::Kokkos::deep_copy(
        ::Kokkos::subview(_ghost_recv_buffer, std::make_pair(std::size_t(0), num_ghosts * stride)),
        ::Kokkos::subview(_ghost_recv_buffer_host, std::make_pair(std::size_t(0), num_ghosts * stride)));
  const auto cloud = _cloud;
  const auto recv = _ghost_recv_buffer;
  ::Kokkos::parallel_for(
      "dem_ghost_forward_unpack", num_ghosts, KOKKOS_LAMBDA(const std::size_t g) {
        for (unsigned int c = 0; c < 3; ++c)
        {
          cloud.x(n + g, c) = recv(g * stride + c);
          cloud.v(n + g, c) = recv(g * stride + 3 + c);
        }
      });
}

void
KokkosParticleCloud::verifyNeighborList()
{
  const auto n = _cloud.size();
  const auto n_total = n + _num_ghosts;
  auto x = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.x);
  auto r = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.r);
  auto gid = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.gid);
  auto offsets =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _neighbor_list.offsets());
  auto pairs = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _neighbor_list.pairs());

  // Pairs with at least one local particle; ghost-ghost pairs are listed but never used
  std::set<std::pair<std::size_t, std::size_t>> listed, expected;
  for (const auto i : make_range(n))
    for (auto k = offsets(i); k < offsets(i + 1); ++k)
      listed.emplace(std::min(i, pairs(k)), std::max(i, pairs(k)));
  for (const auto i : make_range(n))
    for (const auto j : make_range(i + 1, n_total))
    {
      const Point d(x(i, 0) - x(j, 0), x(i, 1) - x(j, 1), x(i, 2) - x(j, 2));
      if (d.norm() < r(i) + r(j) + _skin)
        expected.emplace(i, j);
    }

  for (const auto & [i, j] : expected)
    if (!listed.count({i, j}))
      mooseError("Neighbor list is missing the pair of particles ", gid(i), " and ", gid(j));
  for (const auto & [i, j] : listed)
    if (!expected.count({i, j}))
      mooseError("Neighbor list has the spurious pair of particles ", gid(i), " and ", gid(j));
}

std::size_t
KokkosParticleCloud::migrate()
{
  constexpr auto invalid_proc = libMesh::DofObject::invalid_processor_id;
  constexpr auto exited = DEM::ParticleCloud::exited;
  // Each particle travels as a flat record so no TIMPI type registration is needed; the libMesh
  // element ID it carries is mapped back to a contiguous ID by the receiver
  constexpr auto record_size = DEM::ParticleCloud::record_size;

  // Departing and exited particles are found on device so the common case of none costs one scalar
  const auto cloud = _cloud;
  std::size_t n_sent = 0, n_exited = 0;
  ::Kokkos::parallel_reduce(
      "dem_migrate_count",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sent, std::size_t & gone) {
        sent += cloud.target_rank(i) != invalid_proc;
        gone += cloud.elem(i) == exited;
      },
      n_sent,
      n_exited);
  auto n_removed_global = n_sent + n_exited;
  _communicator.sum(n_removed_global);
  if (!n_removed_global)
    return 0;
  _num_exited_local += n_exited;

  // Prefix sums of the keep and departing masks: a kept particle's slot in the compacted cloud,
  // or a departing particle's slot in the send buffer; exited particles get neither
  const auto slot = _slot;
  ::Kokkos::parallel_scan(
      "dem_migrate_scan_keep",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & n_kept_before, const bool final) {
        const bool keep = cloud.target_rank(i) == invalid_proc && cloud.elem(i) != exited;
        if (final && keep)
          slot(i) = n_kept_before;
        n_kept_before += keep;
      });
  ::Kokkos::parallel_scan(
      "dem_migrate_scan_send",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & n_sent_before, const bool final) {
        const bool send = cloud.target_rank(i) != invalid_proc;
        if (final && send)
          slot(i) = n_sent_before;
        n_sent_before += send;
      });

  // Compact kept particles into the scratch cloud and pack departing ones into the send buffer
  const auto scratch = _scratch;
  const auto libmesh_id = _libmesh_id;
  const auto send_buffer = _send_buffer;
  const auto send_rank = _send_rank;
  ::Kokkos::parallel_for(
      "dem_migrate_pack", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.target_rank(i) != invalid_proc)
        {
          const auto j = slot(i);
          cloud.pack(i, libmesh_id(cloud.elem(i)), &send_buffer(j * record_size));
          send_rank(j) = cloud.target_rank(i);
        }
        else if (cloud.elem(i) != exited)
          scratch.copyFrom(cloud, i, slot(i));
      });
  std::swap(_cloud, _scratch);
  _cloud.n = cloud.size() - n_sent - n_exited;
  // Ghosts were not compacted along; the rebuild after migration replaces them
  _num_ghosts = 0;

  // Only the departing particles cross to the host
  auto send_buffer_host = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{},
      ::Kokkos::subview(send_buffer, std::make_pair(std::size_t(0), n_sent * record_size)));
  auto send_rank_host = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{}, ::Kokkos::subview(send_rank, std::make_pair(std::size_t(0), n_sent)));
  std::map<processor_id_type, std::vector<Real>> to_send;
  for (const auto j : make_range(n_sent))
    to_send[send_rank_host(j)].insert(to_send[send_rank_host(j)].end(),
                                      &send_buffer_host(j * record_size),
                                      &send_buffer_host(j * record_size) + record_size);

  std::vector<Real> received;
  Parallel::push_parallel_vector_data(
      _communicator,
      to_send,
      [&received](processor_id_type, const std::vector<Real> & data)
      { received.insert(received.end(), data.begin(), data.end()); });
  const auto n_received = received.size() / record_size;
  if (!n_received)
    return n_sent;

  // Append the received particles behind the kept ones
  const auto begin = _cloud.size();
  const auto end = begin + n_received;
  if (end > _cloud.capacity())
    reserve(std::max(end, 2 * _cloud.capacity()));

  const auto & mesh = kokkosMesh();
  DEM::ParticleCloud::HostMirror host(n_received);
  for (const auto r : make_range(n_received))
  {
    const Real * record = received.data() + r * record_size;
    host.unpack(record, r);
    const auto * elem = _fe_problem.mesh().getMesh().elem_ptr(record[record_size - 1]);
    mooseAssert(elem->processor_id() == processor_id(), "Received a particle we do not own");
    host.elem(r) = mesh.getContiguousElementID(elem);
  }
  _cloud.copyFrom(host, begin);
  _cloud.n = end;

  return n_sent;
}

void
KokkosParticleCloud::reserve(const std::size_t capacity)
{
  DEM::ParticleCloud cloud;
  cloud.allocate(capacity);
  cloud.copyFrom(_cloud, _cloud.size());
  cloud.n = _cloud.size();
  _cloud = cloud;
  allocateScratch(capacity);
}

void
KokkosParticleCloud::allocateScratch(const std::size_t capacity)
{
  _scratch.allocate(capacity);
  _slot = ::Kokkos::View<std::size_t *>("dem_migrate_slot", capacity);
  // Worst case is every particle departing
  _send_buffer =
      ::Kokkos::View<Real *>("dem_migrate_send", capacity * DEM::ParticleCloud::record_size);
  _send_rank = ::Kokkos::View<processor_id_type *>("dem_migrate_send_rank", capacity);
}

void
KokkosParticleCloud::count()
{
  const auto cloud = _cloud;
  constexpr auto unresolved = DEM::ParticleCloud::unresolved;

  std::size_t n_unresolved = 0;
  unsigned int max_hops_taken = 0;
  Real ke = 0, ke_rot = 0;
  ::Kokkos::parallel_reduce(
      "dem_count_unresolved",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum) {
        sum += cloud.elem(i) == unresolved;
      },
      n_unresolved);
  ::Kokkos::parallel_reduce(
      "dem_count_hops",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, unsigned int & max) {
        max = cloud.hops(i) > max ? cloud.hops(i) : max;
      },
      ::Kokkos::Max<unsigned int>(max_hops_taken));
  ::Kokkos::parallel_reduce(
      "dem_count_ke",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, Real & sum) {
        sum += 0.5 * cloud.m(i) * libMesh::Kokkos::make_vector_ref(cloud.v, i).norm_sq();
      },
      ke);
  ::Kokkos::parallel_reduce(
      "dem_count_ke_rot",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, Real & sum) {
        sum += 0.5 * cloud.inertia(i) * libMesh::Kokkos::make_vector_ref(cloud.omega, i).norm_sq();
      },
      ke_rot);
  // Spin angular momentum; for a sphere the body and space frames coincide
  for (unsigned int c = 0; c < 3; ++c)
    ::Kokkos::parallel_reduce(
        "dem_count_L",
        cloud.size(),
        KOKKOS_LAMBDA(const std::size_t i, Real & sum) {
          sum += cloud.inertia(i) * cloud.omega(i, c);
        },
        _angular_momentum_vector[c]);

  // Local values; finalize() reduces them across ranks every step
  _num_particles = cloud.size();
  _num_exited = _num_exited_local;
  _num_unresolved = n_unresolved;
  _max_hops_taken = max_hops_taken;
  _num_neighbor_pairs = _neighbor_list.numPairs();
  _kinetic_energy = ke;
  _rotational_kinetic_energy = ke_rot;
  _angular_momentum = angularMomentumMagnitude();
}

Real
KokkosParticleCloud::angularMomentumMagnitude() const
{
  return std::sqrt(Utility::pow<2>(_angular_momentum_vector[0]) +
                   Utility::pow<2>(_angular_momentum_vector[1]) +
                   Utility::pow<2>(_angular_momentum_vector[2]));
}

void
KokkosParticleCloud::verify()
{
  const auto num_local = kokkosMesh().getNumLocalElements();
  const auto & ghost_ids = kokkosMesh().getGhostElemIdMapping();
  auto x_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.x);
  auto gid_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.gid);
  auto elem_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.elem);

  for (const auto i : make_range(_cloud.size()))
  {
    const Point p(x_host(i, 0), x_host(i, 1), x_host(i, 2));
    const auto * located = (*_locator)(p);
    const auto elem = elem_host(i);

    if (elem == DEM::ParticleCloud::exited)
      mooseError("Particle ", gid_host(i), " at ", p, " exited but was not removed");
    else if (elem == DEM::ParticleCloud::unresolved)
    {
      // The locator fallback only leaves a particle unresolved beyond the ghost layer
      if (located && (located->processor_id() == processor_id() || ghost_ids.count(located)))
        mooseError(
            "Particle ", gid_host(i), " at ", p, " is unresolved but is in element ", located->id());
    }
    else if (elem >= num_local)
      mooseError("Particle ", gid_host(i), " at ", p, " is still in a ghost element after migration");
    else if (!_cid_to_elem[elem]->contains_point(p))
      mooseError("Particle ",
                 gid_host(i),
                 " at ",
                 p,
                 " was assigned element ",
                 _cid_to_elem[elem]->id(),
                 " which does not contain it");
  }
}

void
KokkosParticleCloud::finalize()
{
  _communicator.sum(_num_particles);
  _communicator.sum(_num_exited);
  _communicator.sum(_num_unresolved);
  _communicator.sum(_num_migrated);
  _communicator.max(_max_hops_taken);
  _communicator.sum(_num_neighbor_pairs);
  _communicator.max(_num_neighbor_list_builds);
  _communicator.max(_num_ghost_forwards);
  _communicator.sum(_kinetic_energy);
  _communicator.sum(_rotational_kinetic_energy);
  for (auto & L : _angular_momentum_vector)
    _communicator.sum(L);
  _angular_momentum = angularMomentumMagnitude();
}
