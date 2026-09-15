#pragma once

#include "KokkosGeneralUserObject.h"
#include "ParticleCloud.h"
#include "NeighborList.h"

#include "libmesh/point_locator_base.h"

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
 * particle that walks out of the mesh or exhausts max_hops is marked lost and keeps integrating.
 *
 * Migration compacts the cloud and packs departing particles on device, so only the departing
 * and arriving particles cross to the host for the TIMPI exchange. GPU-aware MPI (sending the
 * device buffer directly) is deferred until there is a GPU build to measure it on.
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
  std::size_t numLost() const { return _num_lost; }
  std::size_t numMigrated() const { return _num_migrated; }
  unsigned int maxHops() const { return _max_hops_taken; }
  std::size_t numNeighborPairs() const { return _num_neighbor_pairs; }
  std::size_t numNeighborListBuilds() const { return _num_neighbor_list_builds; }
  Real kineticEnergy() const { return _kinetic_energy; }
  Real rotationalKineticEnergy() const { return _rotational_kinetic_energy; }
  Real angularMomentum() const { return _angular_momentum; }
  ///@}

protected:
  /// Zero the force and torque accumulators and apply body forces
  void computeForces();
  /// First half of Velocity Verlet: half-kick the velocities and drift the positions and
  /// orientations by one substep
  void kickDrift();
  /// Second half of Velocity Verlet: half-kick the velocities with the new forces
  void kick();
  /// Re-resolve the containing element of every particle on device, accumulating hops
  void walk();
  /// Rebuild the neighbor list if some particle has moved more than skin / 2 since the last build
  void updateNeighborList();
  /// Send every particle that walked into a ghost element to the rank owning that element and
  /// receive the particles sent here
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
  /// Initial velocity of every particle
  const RealVectorValue & _initial_velocity;
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
  /// Neighbor list of the local particles
  DEM::NeighborList _neighbor_list;
  /// Second cloud of the same capacity that migrate() compacts into before swapping
  DEM::ParticleCloud _scratch;
  ///@{
  /// Work buffers of migrate(), sized to the cloud capacity so no allocation happens per step:
  /// each particle's slot in the compacted cloud or the send buffer, the packed records of the
  /// departing particles, and their destination ranks
  ::Kokkos::View<std::size_t *> _slot;
  ::Kokkos::View<Real *> _send_buffer;
  ::Kokkos::View<processor_id_type *> _send_rank;
  ///@}
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
  std::size_t _num_lost = 0;
  std::size_t _num_migrated = 0;
  unsigned int _max_hops_taken = 0;
  std::size_t _num_neighbor_pairs = 0;
  std::size_t _num_neighbor_list_builds = 0;
  Real _kinetic_energy = 0;
  Real _rotational_kinetic_energy = 0;
  /// Total spin angular momentum vector and its magnitude
  Real _angular_momentum_vector[3] = {0, 0, 0};
  Real _angular_momentum = 0;
};
