#pragma once

#include "KokkosGeneralUserObject.h"
#include "ParticleCloud.h"
#include "NeighborList.h"
#include "LinearSpringDashpot.h"
#include "Hertz.h"
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
 * The normal contact model, linear spring-dashpot or Hertz, is a compile-time policy selected
 * once at setup (plan decision D5): the force kernel is instantiated per model. Walls are fixed
 * planes given by a point and an inward normal each; a particle overlapping a wall gets the same
 * normal contact force as against a sphere at rest of infinite mass.
 *
 * The domain can be periodic in any direction over the extent of the mesh bounding box. A
 * particle that walks out of the mesh is then unresolved rather than exited: it keeps its
 * unwrapped position, just outside its rank's box, until the end of the step, when it is wrapped
 * by the period, placed by the point locator, and migrated to the rank owning the far side, which
 * need not be adjacent. Contact across a periodic face uses the same ghost machinery: every rank
 * sends, to every rank including itself, the particles whose image under each combination of
 * period shifts lies within the cutoff of that rank's box, with the shifted positions, so the
 * neighbor list needs no minimum image convention.
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
  std::size_t numContacts() const { return _num_contacts; }
  Real coordinationNumber() const
  {
    return _num_particles ? 2.0 * _num_contacts / _num_particles : 0;
  }
  Real loadImbalance() const { return _load_imbalance; }
  Real kineticEnergy() const { return _kinetic_energy; }
  Real rotationalKineticEnergy() const { return _rotational_kinetic_energy; }
  Real potentialEnergy() const { return _potential_energy; }
  Real totalEnergy() const
  {
    return _kinetic_energy + _rotational_kinetic_energy + _potential_energy;
  }
  Real linearMomentum() const { return _linear_momentum; }
  Real angularMomentum() const { return _angular_momentum; }
  Real virialPressure() const { return _virial_pressure; }
  ///@}

protected:
  /// Zero the force and torque accumulators and apply body and contact forces to the local
  /// particles listed in a subset, with the selected contact model
  void computeForces(const ::Kokkos::View<std::size_t *> & subset, const std::size_t count);
  /// The above for a given contact model, instantiated per model (plan decision D5)
  template <typename Model>
  void computeForces(const Model & contact,
                     const ::Kokkos::View<std::size_t *> & subset,
                     const std::size_t count);
  /// Count the overlapping pairs and sum their elastic energy and virial with a contact model
  template <typename Model>
  void countContacts(const Model & contact,
                     std::size_t & num_contacts,
                     Real & pe,
                     Real & virial) const;
  /// Build the Hertz model from the input parameters; default-constructed when not selected
  DEM::Hertz makeHertz() const;
  /// Whether the selected contact model applies any force
  bool contactEnabled() const;
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
  /// Locate the unresolved particles with the point locator on host; one found in an element of
  /// another rank, adjacent or not, is flagged for migration there
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
  /// Magnitude of a momentum vector
  static Real magnitude(const Real (&vector)[3]);

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
  /// Normal contact models; the selected one is shared by sphere-sphere and sphere-wall contacts
  enum class ContactModel
  {
    LINEAR_SPRING_DASHPOT,
    HERTZ
  };
  const ContactModel _contact_model;
  const DEM::LinearSpringDashpot _linear;
  const DEM::Hertz _hertz;
  /// Fixed planar walls
  const DEM::AnalyticWalls _walls;
  /// Directions in which the domain is periodic
  const MultiMooseEnum & _periodic;
  ///@{
  /// Lower corner of the mesh bounding box and its extent in each periodic direction, zero in
  /// the others
  Moose::Kokkos::Real3 _domain_min;
  Moose::Kokkos::Real3 _period;
  ///@}
  /// Volume of the mesh bounding box over the mesh dimension, for the virial pressure
  Real _domain_volume = 0;
  /// Number of substeps per MOOSE time step
  const unsigned int _substeps;
  /// Upper bound on face hops per particle per walk
  const unsigned int _max_hops;
  /// Verlet skin distance of the neighbor list
  const Real _skin;
  /// Broad phase of the neighbor list
  const DEM::BroadPhase _broad_phase;
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
  /// libMesh ID of the element of an unresolved particle that the point locator found on another
  /// rank, carried from resolveUnresolved() to migrate() alongside its target rank
  ::Kokkos::View<dof_id_type *> _far_elem;
  /// Number of ghost particles, stored after the _cloud.n local ones
  std::size_t _num_ghosts = 0;
  /// Largest radius over all ranks
  Real _r_max = 0;
  /// Inflated bounding box of every rank's local elements
  std::vector<libMesh::BoundingBox> _rank_boxes;
  /// Ranks whose bounding box comes within the neighbor-list cutoff of ours, possibly through a
  /// period shift, in which case this rank may be its own neighbor
  std::vector<processor_id_type> _neighbor_ranks;
  /// A destination for ghosts: a neighbor rank (an index into _neighbor_ranks) and the period
  /// shift under which our particles come within the cutoff of its box. Sorted by rank, so the
  /// ghosts sent to one rank are contiguous in the send index
  struct GhostTarget
  {
    std::size_t neighbor;
    Moose::Kokkos::Real3 shift;
  };
  std::vector<GhostTarget> _ghost_targets;
  ///@{
  /// Ghost forwarding state fixed between rebuilds: the local particles ghosted to each neighbor
  /// rank in _neighbor_ranks order, concatenated, with their counts; the source rank, first ghost
  /// slot, and count of each received block; the device send and receive buffers with host
  /// mirrors for MPI libraries that cannot read device memory; the pending requests; and the
  /// message tag
  ::Kokkos::View<std::size_t *> _ghost_send_index;
  /// Period shift applied to the positions of each ghosted particle in _ghost_send_index
  ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight> _ghost_send_shift;
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
  /// Overlapping pairs, each counted once across ranks
  std::size_t _num_contacts = 0;
  /// Largest local particle count over the mean, 1 when balanced
  Real _load_imbalance = 0;
  Real _kinetic_energy = 0;
  Real _rotational_kinetic_energy = 0;
  /// Elastic energy of the overlapping pairs and wall contacts
  Real _potential_energy = 0;
  /// Pair virial sum over the overlapping pairs, r_ij . f_ij, each counted once
  Real _virial = 0;
  Real _virial_pressure = 0;
  ///@{
  /// Total linear and spin angular momentum vectors and their magnitudes
  Real _linear_momentum_vector[3] = {0, 0, 0};
  Real _linear_momentum = 0;
  Real _angular_momentum_vector[3] = {0, 0, 0};
  Real _angular_momentum = 0;
  ///@}
};
