#pragma once

#include "KokkosGeneralUserObject.h"
#include "ParticleCloud.h"
#include "NeighborList.h"
#include "LinearSpringDashpot.h"
#include "AnalyticWalls.h"

#include "libmesh/point_locator_base.h"
#include "libmesh/bounding_box.h"
#include "libmesh/parallel.h"

#include <mpi.h>

/**
 * A cloud of spherical particles integrated on device with Velocity Verlet over a number of
 * substeps per MOOSE time step (plan decision D4), tracking the local element containing each
 * particle by walking across element faces (decision D3, layer L1).
 *
 * The walk uses the cached side centroids and outward normals of Moose::Kokkos::Mesh, so it is
 * exact only for elements with planar faces (TRI3, QUAD4, TET4, HEX8, ...). A particle that walks
 * into a ghost element keeps integrating with a stale element until the end of the MOOSE step,
 * when it is sent to the rank owning that element, which resumes the walk; this repeats until no
 * particle is left in a ghost element, so a particle may cross several partitions in one step. A
 * particle that walks out of the mesh has exited: it is inert for the rest of the step and removed
 * at the end. A particle that exhausts max_hops is unresolved: it keeps integrating and is located
 * again by a point locator at the end of the step.
 *
 * Contact across partitions uses ghost particles: copies of the neighboring ranks' particles
 * within the neighbor-list cutoff of this rank's bounding box, appended after the local particles.
 * The ghost list is rebuilt with the neighbor list and ghost positions and velocities are forwarded
 * every substep with nonblocking MPI directly from the device buffers when the host can access
 * device memory or the MPI library is GPU-aware, and through host mirrors otherwise. Forces on
 * particles without ghost neighbors are computed while the forward is in flight. Every rank
 * computes the forces on its own particles from the ghost copies, so no force communication is
 * needed.
 *
 * Walls are fixed planes given by a point and an inward normal each; a particle overlapping a
 * wall gets the same normal contact force as against a sphere at rest of infinite mass.
 *
 * Migration packs on device and stages through host buffers for the TIMPI exchange.
 */
class KokkosParticleCloud : public Moose::Kokkos::GeneralUserObject
{
public:
  static InputParameters validParams();

  KokkosParticleCloud(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override {}
  virtual void compute() override;
  virtual void finalize() override;

  /// Device particle state
  const DEM::ParticleCloud & cloud() const { return _cloud; }

  ///@{
  /// Counters, global across ranks after finalize()
  std::size_t numParticles() const { return _num_particles; }
  std::size_t numExited() const { return _num_exited; }
  std::size_t numUnresolved() const { return _num_unresolved; }
  std::size_t numMigrated() const { return _num_migrated; }
  unsigned int maxHops() const { return _max_hops_taken; }
  std::size_t numNeighborPairs() const { return _num_neighbor_pairs; }
  std::size_t numNeighborListBuilds() const { return _num_neighbor_list_builds; }
  std::size_t numGhosts() const { return _num_ghosts; }
  std::size_t numGhostForwards() const { return _num_ghost_forwards; }
  Real kineticEnergy() const { return _kinetic_energy; }
  Real rotationalKineticEnergy() const { return _rotational_kinetic_energy; }
  Real angularMomentum() const { return _angular_momentum; }
  ///@}

protected:
  /// Zero the force and torque accumulators and apply body and contact forces to the local
  /// particles listed in a subset
  void computeForces(const ::Kokkos::View<std::size_t *> & subset, const std::size_t count);
  /// Compute all forces: interior particles first, then the boundary ones once the ghost update
  /// has landed
  void computeForces();
  /// First half of Velocity Verlet: half-kick the velocities and drift the positions and
  /// orientations by one substep
  void kickDrift();
  /// Second half of Velocity Verlet: half-kick the velocities with the new forces
  void kick();
  /// Re-resolve the containing element of every particle on device, accumulating hops
  void walk();
  /// Locate the unresolved particles with the point locator on host
  void resolveUnresolved();
  /// Rebuild the ghost and neighbor lists if, on any rank, some particle has moved more than
  /// skin / 2 since the last build or the local particle count changed
  /// @returns Whether the lists were rebuilt
  bool updateNeighborList();
  /// Replace the ghost particles with the neighboring ranks' particles within the neighbor-list
  /// cutoff of this rank's bounding box, and record what to forward every substep
  void exchangeGhosts();
  /// Start forwarding the current positions and velocities of the ghosted particles
  void startGhostUpdate();
  /// Wait for the ghost update started by startGhostUpdate() and unpack it; no-op otherwise
  void finishGhostUpdate();
  /// Remove the exited particles, send every particle that walked into a ghost element to the
  /// rank owning that element, and receive the particles sent here
  /// @returns The number of particles sent from this rank
  std::size_t migrate();
  /// Reallocate both clouds and the migration work buffers so they can hold at least the given
  /// number of particles, keeping the live particles of _cloud
  void reserve(const std::size_t capacity);
  /// Allocate the scratch cloud and migration work buffers for a capacity
  void allocateScratch(const std::size_t capacity);
  /// Update the local counters and energies from the device state
  void count();
  /// Check the device element assignment against libMesh's PointLocator on host; errors on mismatch
  void verify();
  /// Check the neighbor list against a brute-force pair search on host; errors on mismatch
  void verifyNeighborList();
  /// Magnitude of _angular_momentum_vector
  Real angularMomentumMagnitude() const;

  /// Initial particle positions; a rank keeps only the particles in its local elements
  const std::vector<Point> & _initial_positions;
  /// Initial velocity of every particle, used unless initial_velocities is given
  const RealVectorValue & _initial_velocity;
  /// Optional per-particle initial velocities
  const std::vector<Point> _initial_velocities;
  /// Initial body-frame angular velocity of every particle
  const RealVectorValue & _initial_angular_velocity;
  /// Particle radius, used unless initial_radii is given
  const Real _radius;
  /// Optional per-particle radii
  const std::vector<Real> _initial_radii;
  /// Particle density
  const Real _density;
  /// Gravitational acceleration
  const RealVectorValue & _gravity;
  /// Normal contact model, shared by sphere-sphere and sphere-wall contacts
  const DEM::LinearSpringDashpot _contact;
  /// Fixed planar walls
  const DEM::AnalyticWalls _walls;
  /// Number of substeps per MOOSE time step
  const unsigned int _substeps;
  /// Upper bound on face hops per particle per walk
  const unsigned int _max_hops;
  /// Verlet skin distance of the neighbor list
  const Real _skin;
  /// Whether to run verify() every step
  const bool _verify;
  /// Whether to run verifyNeighborList() after every build
  const bool _verify_neighbor_list;

  /// Device particle state
  DEM::ParticleCloud _cloud;
  /// Neighbor list of the local and ghost particles
  DEM::NeighborList _neighbor_list;
  /// Second cloud of the same capacity that migrate() compacts into before swapping
  DEM::ParticleCloud _scratch;
  ///@{
  /// Work buffers of migrate(), sized to the cloud capacity so no allocation happens per step:
  /// each particle's slot in the compacted cloud or the send buffer (also used to flag ghost
  /// candidates), the packed records of the departing particles, and their destination ranks
  ::Kokkos::View<std::size_t *> _slot;
  ::Kokkos::View<Real *> _send_buffer;
  ::Kokkos::View<processor_id_type *> _send_rank;
  ///@}
  /// Number of ghost particles, stored after the _cloud.n local ones
  std::size_t _num_ghosts = 0;
  /// Largest radius over all ranks
  Real _r_max = 0;
  /// Inflated bounding box of every rank's local elements
  std::vector<libMesh::BoundingBox> _rank_boxes;
  /// Ranks whose bounding box comes within the neighbor-list cutoff of ours
  std::vector<processor_id_type> _neighbor_ranks;
  ///@{
  /// Ghost forwarding state fixed between rebuilds: the local particles ghosted to each neighbor
  /// rank in _neighbor_ranks order, concatenated, with their counts; the source rank, first ghost
  /// slot, and count of each received block; the device send and receive buffers with host
  /// mirrors for MPI libraries that cannot read device memory; the pending requests; and the
  /// message tag
  ::Kokkos::View<std::size_t *> _ghost_send_index;
  std::vector<std::size_t> _ghost_send_counts;
  struct GhostBlock
  {
    processor_id_type source;
    std::size_t begin;
    std::size_t count;
  };
  std::vector<GhostBlock> _ghost_recv;
  ::Kokkos::View<Real *> _ghost_buffer;
  ::Kokkos::View<Real *> _ghost_recv_buffer;
  ::Kokkos::View<Real *>::HostMirror _ghost_buffer_host;
  ::Kokkos::View<Real *>::HostMirror _ghost_recv_buffer_host;
  std::vector<MPI_Request> _ghost_requests;
  const libMesh::Parallel::MessageTag _ghost_tag;
  ///@}
  /// Whether the ghost forward may pass device buffers to MPI
  static constexpr bool _device_buffers_to_mpi =
      ::Kokkos::SpaceAccessibility<::Kokkos::HostSpace,
                                   ::Kokkos::DefaultExecutionSpace::memory_space>::accessible;
  /// Whether to pass device buffers to MPI on a GPU build (requires a GPU-aware MPI library)
  const bool _gpu_aware_mpi;
  /// Owning rank of every local and one-layer ghost element, indexed by contiguous element ID
  ::Kokkos::View<processor_id_type *> _owner;
  /// libMesh element ID of every local and one-layer ghost element, indexed by contiguous element
  /// ID; the rank-independent ID a migrating particle carries
  ::Kokkos::View<dof_id_type *> _libmesh_id;
  /// Reverse of Moose::Kokkos::Mesh::getContiguousElementID for local and ghost elements
  std::vector<const Elem *> _cid_to_elem;
  /// Point locator for initial placement and verification
  std::unique_ptr<libMesh::PointLocatorBase> _locator;

  std::size_t _num_particles = 0;
  /// Particles removed after walking out of the mesh, cumulative on this rank and globally
  std::size_t _num_exited_local = 0;
  std::size_t _num_exited = 0;
  std::size_t _num_unresolved = 0;
  std::size_t _num_migrated = 0;
  unsigned int _max_hops_taken = 0;
  std::size_t _num_neighbor_pairs = 0;
  std::size_t _num_neighbor_list_builds = 0;
  std::size_t _num_ghost_forwards = 0;
  Real _kinetic_energy = 0;
  Real _rotational_kinetic_energy = 0;
  /// Total spin angular momentum vector and its magnitude
  Real _angular_momentum_vector[3] = {0, 0, 0};
  Real _angular_momentum = 0;
};
