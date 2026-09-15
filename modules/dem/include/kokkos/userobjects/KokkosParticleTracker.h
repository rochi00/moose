#pragma once

#include "KokkosGeneralUserObject.h"
#include "ParticleCloud.h"

#include "libmesh/point_locator_base.h"

/**
 * Advects a cloud of particles with a uniform velocity and tracks the local element containing
 * each particle by walking across element faces on device (plan decision D3, layer L1).
 *
 * The walk uses the cached side centroids and outward normals of Moose::Kokkos::Mesh, so it is
 * exact only for elements with planar faces (TRI3, QUAD4, TET4, HEX8, ...). A particle that walks
 * into a ghost element is sent to the rank owning that element, which resumes the walk; this
 * repeats until no particle is left in a ghost element, so a particle may cross several partitions
 * in one step. A particle that walks out of the mesh or exhausts max_hops is marked lost.
 *
 * Migration compacts the cloud and packs departing particles on device, so only the departing
 * and arriving particles cross to the host for the TIMPI exchange. GPU-aware MPI (sending the
 * device buffer directly) is deferred until there is a GPU build to measure it on.
 */
class KokkosParticleTracker : public Moose::Kokkos::GeneralUserObject
{
public:
  static InputParameters validParams();

  KokkosParticleTracker(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override {}
  virtual void compute() override;
  virtual void finalize() override;

  ///@{
  /// Counters, global across ranks after finalize()
  std::size_t numParticles() const { return _num_particles; }
  std::size_t numLost() const { return _num_lost; }
  std::size_t numMigrated() const { return _num_migrated; }
  unsigned int maxHops() const { return _max_hops_taken; }
  ///@}

protected:
  /// Advance positions by _dt on device
  void advance();
  /// Re-resolve the containing element of every particle on device, accumulating hops
  void walk();
  /// Send every particle that walked into a ghost element to the rank owning that element and
  /// receive the particles sent here
  /// @returns The number of particles sent from this rank
  std::size_t migrate();
  /// Reallocate both clouds and the migration work buffers so they can hold at least the given
  /// number of particles, keeping the live particles of _cloud
  void reserve(const std::size_t capacity);
  /// Allocate the scratch cloud and migration work buffers for a capacity
  void allocateScratch(const std::size_t capacity);
  /// Update the local counters from the device state
  void count();
  /// Check the device element assignment against libMesh's PointLocator on host; errors on mismatch
  void verify();

  /// Initial particle positions; a rank keeps only the particles in its local elements
  const std::vector<Point> & _initial_positions;
  /// Uniform particle velocity
  const RealVectorValue & _velocity;
  /// Upper bound on face hops per particle per step
  const unsigned int _max_hops;
  /// Whether to run verify() every step
  const bool _verify;

  /// Device particle state
  DEM::ParticleCloud _cloud;
  /// Second cloud of the same capacity that migrate() compacts into before swapping
  DEM::ParticleCloud _scratch;
  /// Number of Reals a migrating particle is packed into: gid, x, v, hops, libMesh element ID
  static constexpr std::size_t _record_size = 9;
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
};
