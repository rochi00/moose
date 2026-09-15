#include "KokkosParticleTracker.h"

#include "MooseMesh.h"

#include "libmesh/elem.h"
#include "timpi/parallel_sync.h"

registerMooseObject("DemApp", KokkosParticleTracker);

InputParameters
KokkosParticleTracker::validParams()
{
  InputParameters params = Moose::Kokkos::GeneralUserObject::validParams();
  params.addClassDescription("Advects particles with a uniform velocity and tracks the element "
                             "containing each particle by walking across element faces on device, "
                             "migrating particles between ranks as they cross partitions.");
  params.addRequiredParam<std::vector<Point>>("initial_positions",
                                              "Initial positions of all particles");
  params.addParam<RealVectorValue>(
      "velocity", RealVectorValue(0, 0, 0), "Uniform velocity applied to every particle");
  // A particle should cross only a few cells per step; a much larger walk indicates a mesh the
  // face-plane test cannot handle, so bound it rather than loop for the mesh size
  params.addParam<unsigned int>(
      "max_hops", 32, "Maximum number of face-neighbor hops per particle per step");
  params.addParam<bool>("verify",
                        false,
                        "Check the device element assignment against libMesh's PointLocator "
                        "every step and error on any mismatch");
  return params;
}

KokkosParticleTracker::KokkosParticleTracker(const InputParameters & parameters)
  : Moose::Kokkos::GeneralUserObject(parameters),
    _initial_positions(getParam<std::vector<Point>>("initial_positions")),
    _velocity(getParam<RealVectorValue>("velocity")),
    _max_hops(getParam<unsigned int>("max_hops")),
    _verify(getParam<bool>("verify"))
{
  // Side centroids and normals are only cached on request
  _fe_problem.mesh().getKokkosMesh()->setNeedsElementSideGeometry();
}

void
KokkosParticleTracker::initialSetup()
{
  const auto & mesh = kokkosMesh();
  const auto num_total = mesh.getNumLocalAndPossiblyOneNeighborLayerGhostElements();

  // Reverse element ID map and owner ranks for local and one-layer ghost elements
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

  std::vector<int64_t> gids;
  std::vector<Point> positions;
  std::vector<ContiguousElementID> elems;
  for (const auto i : index_range(_initial_positions))
  {
    const auto * elem = (*_locator)(_initial_positions[i]);
    if (!elem || elem->processor_id() != processor_id())
      continue;
    gids.push_back(i);
    positions.push_back(_initial_positions[i]);
    elems.push_back(mesh.getContiguousElementID(elem));
  }

  const auto n = gids.size();
  _cloud.allocate(n);
  _cloud.x = libMesh::Kokkos::upload_vector_storage<DEM::ParticleCloud::StoragePolicy>(positions,
                                                                                       "dem_x");
  _cloud.v = libMesh::Kokkos::upload_vector_storage<DEM::ParticleCloud::StoragePolicy>(
      std::vector<RealVectorValue>(n, _velocity), "dem_v");

  auto gid_host = ::Kokkos::create_mirror_view(_cloud.gid);
  auto elem_host = ::Kokkos::create_mirror_view(_cloud.elem);
  for (const auto i : make_range(n))
  {
    gid_host(i) = gids[i];
    elem_host(i) = elems[i];
  }
  ::Kokkos::deep_copy(_cloud.gid, gid_host);
  ::Kokkos::deep_copy(_cloud.elem, elem_host);
  ::Kokkos::deep_copy(_cloud.target_rank, libMesh::DofObject::invalid_processor_id);
  ::Kokkos::deep_copy(_cloud.hops, 0u);
  _cloud.n = n;
  allocateScratch(n);
}

void
KokkosParticleTracker::compute()
{
  advance();

  // Walk, hand particles that reached a ghost element to their owner, and let the owner continue
  // the walk until every particle has settled. Each round moves a particle at least one partition,
  // so the number of rounds is bounded by the number of ranks
  _num_migrated = 0;
  for (const auto round : make_range(n_processors()))
  {
    walk();

    auto n_sent = migrate();
    _num_migrated += n_sent;
    _communicator.sum(n_sent);
    if (!n_sent)
      break;

    if (round + 1 == n_processors())
      mooseError("Particles were still crossing partitions after ", n_processors(), " rounds");
  }

  count();
  if (_verify)
    verify();
}

void
KokkosParticleTracker::advance()
{
  const auto cloud = _cloud;
  const Real dt = _dt;

  ::Kokkos::parallel_for(
      "dem_advance", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        auto x = libMesh::Kokkos::make_vector_ref(cloud.x, i);
        x.add_scaled(libMesh::Kokkos::make_vector_ref(cloud.v, i), dt);
        cloud.hops(i) = 0;
      });
}

void
KokkosParticleTracker::walk()
{
  // Shallow copies captured by value into the kernel
  const Moose::Kokkos::Mesh mesh = kokkosMesh();
  const auto num_local = mesh.getNumLocalElements();
  const auto cloud = _cloud;
  const auto owner = _owner;
  const auto max_hops = _max_hops;
  constexpr auto invalid_elem = libMesh::DofObject::invalid_id;
  constexpr auto invalid_proc = libMesh::DofObject::invalid_processor_id;

  ::Kokkos::parallel_for(
      "dem_walk", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        auto elem = cloud.elem(i);
        cloud.target_rank(i) = invalid_proc;
        if (elem == invalid_elem)
          return;

        const auto x = libMesh::Kokkos::make_vector_ref(cloud.x, i);
        const Moose::Kokkos::Real3 xp(x(0), x(1), x(2));
        auto hops = cloud.hops(i);
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
            elem = invalid_elem;
            break;
          }
          elem = mesh.getNeighbor(elem, exit_side);
          ++hops;
          if (elem == invalid_elem)
            break;
        }
        cloud.elem(i) = elem;
        cloud.hops(i) = hops;
      });
}

std::size_t
KokkosParticleTracker::migrate()
{
  constexpr auto invalid_proc = libMesh::DofObject::invalid_processor_id;
  // Each particle travels as a flat record so no TIMPI type registration is needed; the libMesh
  // element ID it carries is mapped back to a contiguous ID by the receiver
  constexpr auto record_size = _record_size;

  // Departing particles are found on device so the common case of none costs one scalar
  const auto cloud = _cloud;
  std::size_t n_sent = 0;
  ::Kokkos::parallel_reduce(
      "dem_migrate_count",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum) {
        sum += cloud.target_rank(i) != invalid_proc;
      },
      n_sent);
  auto n_sent_global = n_sent;
  _communicator.sum(n_sent_global);
  if (!n_sent_global)
    return 0;

  // Prefix sum of the keep mask: a kept particle's slot in the compacted cloud, or a departing
  // particle's slot in the send buffer
  const auto slot = _slot;
  ::Kokkos::parallel_scan(
      "dem_migrate_scan",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & n_kept_before, const bool final) {
        const bool keep = cloud.target_rank(i) == invalid_proc;
        if (final)
          slot(i) = keep ? n_kept_before : i - n_kept_before;
        n_kept_before += keep;
      });

  // Compact kept particles into the scratch cloud and pack departing ones into the send buffer
  const auto scratch = _scratch;
  const auto libmesh_id = _libmesh_id;
  const auto send_buffer = _send_buffer;
  const auto send_rank = _send_rank;
  ::Kokkos::parallel_for(
      "dem_migrate_pack", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        const auto j = slot(i);
        if (cloud.target_rank(i) == invalid_proc)
        {
          scratch.copyFrom(cloud, i, j);
          return;
        }
        Real * record = &send_buffer(j * record_size);
        record[0] = cloud.gid(i);
        for (unsigned int c = 0; c < 3; ++c)
        {
          record[1 + c] = cloud.x(i, c);
          record[4 + c] = cloud.v(i, c);
        }
        record[7] = cloud.hops(i);
        record[8] = libmesh_id(cloud.elem(i));
        send_rank(j) = cloud.target_rank(i);
      });
  std::swap(_cloud, _scratch);
  _cloud.n = cloud.size() - n_sent;

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
  auto gid_host = ::Kokkos::View<int64_t *, ::Kokkos::HostSpace>("dem_recv_gid", n_received);
  auto x_host = ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight, ::Kokkos::HostSpace>(
      "dem_recv_x", n_received);
  auto v_host = ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight, ::Kokkos::HostSpace>(
      "dem_recv_v", n_received);
  auto elem_host =
      ::Kokkos::View<ContiguousElementID *, ::Kokkos::HostSpace>("dem_recv_elem", n_received);
  auto hops_host = ::Kokkos::View<unsigned int *, ::Kokkos::HostSpace>("dem_recv_hops", n_received);
  for (const auto r : make_range(n_received))
  {
    const Real * record = received.data() + r * record_size;
    gid_host(r) = record[0];
    for (const auto c : make_range(3))
    {
      x_host(r, c) = record[1 + c];
      v_host(r, c) = record[4 + c];
    }
    hops_host(r) = record[7];
    const auto * elem = _fe_problem.mesh().getMesh().elem_ptr(record[8]);
    mooseAssert(elem->processor_id() == processor_id(), "Received a particle we do not own");
    elem_host(r) = mesh.getContiguousElementID(elem);
  }
  const auto range = std::make_pair(begin, end);
  ::Kokkos::deep_copy(::Kokkos::subview(_cloud.gid, range), gid_host);
  ::Kokkos::deep_copy(::Kokkos::subview(_cloud.x, range, ::Kokkos::ALL), x_host);
  ::Kokkos::deep_copy(::Kokkos::subview(_cloud.v, range, ::Kokkos::ALL), v_host);
  ::Kokkos::deep_copy(::Kokkos::subview(_cloud.elem, range), elem_host);
  ::Kokkos::deep_copy(::Kokkos::subview(_cloud.hops, range), hops_host);
  ::Kokkos::deep_copy(::Kokkos::subview(_cloud.target_rank, range), invalid_proc);
  _cloud.n = end;

  return n_sent;
}

void
KokkosParticleTracker::reserve(const std::size_t capacity)
{
  DEM::ParticleCloud cloud;
  cloud.allocate(capacity);
  const auto range = std::make_pair(std::size_t(0), _cloud.size());
  ::Kokkos::deep_copy(::Kokkos::subview(cloud.gid, range), ::Kokkos::subview(_cloud.gid, range));
  ::Kokkos::deep_copy(::Kokkos::subview(cloud.x, range, ::Kokkos::ALL),
                      ::Kokkos::subview(_cloud.x, range, ::Kokkos::ALL));
  ::Kokkos::deep_copy(::Kokkos::subview(cloud.v, range, ::Kokkos::ALL),
                      ::Kokkos::subview(_cloud.v, range, ::Kokkos::ALL));
  ::Kokkos::deep_copy(::Kokkos::subview(cloud.elem, range), ::Kokkos::subview(_cloud.elem, range));
  ::Kokkos::deep_copy(::Kokkos::subview(cloud.target_rank, range),
                      ::Kokkos::subview(_cloud.target_rank, range));
  ::Kokkos::deep_copy(::Kokkos::subview(cloud.hops, range), ::Kokkos::subview(_cloud.hops, range));
  cloud.n = _cloud.size();
  _cloud = cloud;
  allocateScratch(capacity);
}

void
KokkosParticleTracker::allocateScratch(const std::size_t capacity)
{
  _scratch.allocate(capacity);
  _slot = ::Kokkos::View<std::size_t *>("dem_migrate_slot", capacity);
  // Worst case is every particle departing
  _send_buffer = ::Kokkos::View<Real *>("dem_migrate_send", capacity * _record_size);
  _send_rank = ::Kokkos::View<processor_id_type *>("dem_migrate_send_rank", capacity);
}

void
KokkosParticleTracker::count()
{
  const auto cloud = _cloud;
  constexpr auto invalid_elem = libMesh::DofObject::invalid_id;

  std::size_t lost = 0;
  unsigned int max_hops_taken = 0;
  ::Kokkos::parallel_reduce(
      "dem_count_lost",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum) {
        sum += cloud.elem(i) == invalid_elem;
      },
      lost);
  ::Kokkos::parallel_reduce(
      "dem_count_hops",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, unsigned int & max) {
        max = cloud.hops(i) > max ? cloud.hops(i) : max;
      },
      ::Kokkos::Max<unsigned int>(max_hops_taken));

  // Local counts; finalize() reduces them across ranks every step
  _num_particles = cloud.size();
  _num_lost = lost;
  _max_hops_taken = max_hops_taken;
}

void
KokkosParticleTracker::verify()
{
  const auto num_local = kokkosMesh().getNumLocalElements();
  auto x_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.x);
  auto gid_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.gid);
  auto elem_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.elem);

  for (const auto i : make_range(_cloud.size()))
  {
    const Point p(x_host(i, 0), x_host(i, 1), x_host(i, 2));
    const auto * located = (*_locator)(p);
    const auto elem = elem_host(i);

    if (elem == libMesh::DofObject::invalid_id)
    {
      // Lost: the point must really be outside the mesh
      if (located)
        mooseError("Particle ", gid_host(i), " at ", p, " was lost but is in element ", located->id());
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
KokkosParticleTracker::finalize()
{
  _communicator.sum(_num_particles);
  _communicator.sum(_num_lost);
  _communicator.sum(_num_migrated);
  _communicator.max(_max_hops_taken);
}
