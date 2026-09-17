#pragma once

#include "KokkosMesh.h"

#include "libmesh/kokkos_storage_policy.h"
#include "libmesh/kokkos_storage.h"
#include "libmesh/kokkos_vector_ops.h"

namespace DEM
{

/**
 * Device-resident particle state (layer L0 of the DEM module plan).
 *
 * Vector quantities use the libMesh Kokkos storage policy so kernels can operate on entries
 * through libMesh::Kokkos::vector_ref / Point without a third vector type. LayoutRight keeps each
 * particle's components contiguous, which matches Moose::Kokkos::Real3 and LAMMPS-KOKKOS; the
 * policy is a single alias so it can be benchmarked against LayoutLeft later (plan Section 4.4).
 *
 * Element IDs are the Kokkos-MOOSE contiguous local element IDs (Moose::Kokkos::Mesh), not libMesh
 * IDs, and are therefore only meaningful on the rank that owns the particle.
 */
struct ParticleCloud
{
  using StoragePolicy = libMesh::Kokkos::layout_right_storage_policy;
  using VectorView = StoragePolicy::vector_view;

  /// Global particle ID, unique for the lifetime of the run (plan decision D2)
  ::Kokkos::View<int64_t *> gid;
  /// Position
  VectorView x;
  /// Velocity
  VectorView v;
  /// Angular velocity in the body frame
  VectorView omega;
  /// Orientation quaternion (w, x, y, z) rotating the body frame into the space frame
  ::Kokkos::View<Real * [4], ::Kokkos::LayoutRight> q;
  /// Radius
  ::Kokkos::View<Real *> r;
  /// Mass
  ::Kokkos::View<Real *> m;
  /// Moment of inertia (scalar; spheres only)
  ::Kokkos::View<Real *> inertia;
  /// Force accumulator, zeroed every substep; between steps it holds the force of the last
  /// substep, which the next step's first half-kick uses
  VectorView f;
  /// Torque accumulator in the body frame, zeroed every substep
  VectorView tau;
  /// Force and torque recomputed at the end of the step from the settled positions and the
  /// full-step velocities, for output; the half-kick uses f and tau, computed with the
  /// half-step velocities like every substep's, so the trajectory does not depend on how the
  /// substeps are grouped into steps
  VectorView f_out;
  VectorView tau_out;
  /// Contiguous local element containing the particle, or one of the sentinels below
  ::Kokkos::View<ContiguousElementID *> elem;
  /// elem value of a particle whose element could not be resolved by the face walk; retried by a
  /// point locator at the end of the step
  static constexpr ContiguousElementID unresolved = libMesh::DofObject::invalid_id;
  /// elem value of a particle that walked out of the mesh; inert until removed at the end of the
  /// step
  static constexpr ContiguousElementID exited = libMesh::DofObject::invalid_id - 1;
  /// Rank that should receive the particle when it has walked into a ghost element;
  /// libMesh::DofObject::invalid_processor_id when the particle is local or lost
  ::Kokkos::View<processor_id_type *> target_rank;
  /// Number of face-neighbor hops the tracker took for this particle in the current step
  ::Kokkos::View<unsigned int *> hops;
  /// Nonzero for a frozen particle: never integrated, an immovable obstacle to the others
  ::Kokkos::View<unsigned char *> frozen;

  /// Number of live particles; entries [n, capacity()) are unused
  std::size_t n = 0;

  /// Number of live particles
  std::size_t size() const { return n; }
  /// Number of particles the views can hold
  std::size_t capacity() const { return gid.extent(0); }

  /// Allocate all views for a given capacity with no live particles (contents uninitialized)
  void allocate(const std::size_t capacity)
  {
    gid = ::Kokkos::View<int64_t *>("dem_gid", capacity);
    x = libMesh::Kokkos::make_vector_storage<StoragePolicy>("dem_x", capacity);
    v = libMesh::Kokkos::make_vector_storage<StoragePolicy>("dem_v", capacity);
    omega = libMesh::Kokkos::make_vector_storage<StoragePolicy>("dem_omega", capacity);
    q = ::Kokkos::View<Real * [4], ::Kokkos::LayoutRight>("dem_q", capacity);
    r = ::Kokkos::View<Real *>("dem_r", capacity);
    m = ::Kokkos::View<Real *>("dem_m", capacity);
    inertia = ::Kokkos::View<Real *>("dem_inertia", capacity);
    f = libMesh::Kokkos::make_vector_storage<StoragePolicy>("dem_f", capacity);
    tau = libMesh::Kokkos::make_vector_storage<StoragePolicy>("dem_tau", capacity);
    f_out = libMesh::Kokkos::make_vector_storage<StoragePolicy>("dem_f_out", capacity);
    tau_out = libMesh::Kokkos::make_vector_storage<StoragePolicy>("dem_tau_out", capacity);
    elem = ::Kokkos::View<ContiguousElementID *>("dem_elem", capacity);
    target_rank = ::Kokkos::View<processor_id_type *>("dem_target_rank", capacity);
    hops = ::Kokkos::View<unsigned int *>("dem_hops", capacity);
    frozen = ::Kokkos::View<unsigned char *>("dem_frozen", capacity);
    n = 0;
  }

  /// Copy particle i of another cloud into slot j of this one
  KOKKOS_INLINE_FUNCTION void copyFrom(const ParticleCloud & other,
                                       const std::size_t i,
                                       const std::size_t j) const
  {
    gid(j) = other.gid(i);
    for (unsigned int c = 0; c < 3; ++c)
    {
      x(j, c) = other.x(i, c);
      v(j, c) = other.v(i, c);
      omega(j, c) = other.omega(i, c);
      f(j, c) = other.f(i, c);
      tau(j, c) = other.tau(i, c);
    }
    for (unsigned int c = 0; c < 4; ++c)
      q(j, c) = other.q(i, c);
    r(j) = other.r(i);
    m(j) = other.m(i);
    inertia(j) = other.inertia(i);
    elem(j) = other.elem(i);
    target_rank(j) = other.target_rank(i);
    hops(j) = other.hops(i);
    frozen(j) = other.frozen(i);
  }

  /// Number of Reals a particle is packed into for migration: gid, x, v, omega, q, r, m, inertia,
  /// hops, f, tau, the frozen flag, and the libMesh ID of its element
  static constexpr std::size_t record_size = 26;

  /// Pack particle i into a migration record, with its element given as a libMesh ID
  KOKKOS_INLINE_FUNCTION void pack(const std::size_t i,
                                   const dof_id_type libmesh_elem_id,
                                   Real * const record) const
  {
    record[0] = gid(i);
    for (unsigned int c = 0; c < 3; ++c)
    {
      record[1 + c] = x(i, c);
      record[4 + c] = v(i, c);
      record[7 + c] = omega(i, c);
    }
    for (unsigned int c = 0; c < 4; ++c)
      record[10 + c] = q(i, c);
    record[14] = r(i);
    record[15] = m(i);
    record[16] = inertia(i);
    record[17] = hops(i);
    for (unsigned int c = 0; c < 3; ++c)
    {
      record[18 + c] = f(i, c);
      record[21 + c] = tau(i, c);
    }
    record[24] = frozen(i);
    record[25] = libmesh_elem_id;
  }

  /// Host-side copy of the particle views, used to stage received particles
  struct HostMirror
  {
    ::Kokkos::View<int64_t *>::HostMirror gid;
    VectorView::HostMirror x, v, omega, f, tau;
    ::Kokkos::View<Real * [4], ::Kokkos::LayoutRight>::HostMirror q;
    ::Kokkos::View<Real *>::HostMirror r, m, inertia;
    ::Kokkos::View<ContiguousElementID *>::HostMirror elem;
    ::Kokkos::View<processor_id_type *>::HostMirror target_rank;
    ::Kokkos::View<unsigned int *>::HostMirror hops;
    ::Kokkos::View<unsigned char *>::HostMirror frozen;

    HostMirror(const std::size_t n)
      : gid("dem_host_gid", n),
        x("dem_host_x", n),
        v("dem_host_v", n),
        omega("dem_host_omega", n),
        f("dem_host_f", n),
        tau("dem_host_tau", n),
        q("dem_host_q", n),
        r("dem_host_r", n),
        m("dem_host_m", n),
        inertia("dem_host_inertia", n),
        elem("dem_host_elem", n),
        target_rank("dem_host_target_rank", n),
        hops("dem_host_hops", n),
        frozen("dem_host_frozen", n)
    {
    }

    /// Unpack a migration record into slot j; the element is left to the caller
    void unpack(const Real * const record, const std::size_t j) const
    {
      gid(j) = record[0];
      for (unsigned int c = 0; c < 3; ++c)
      {
        x(j, c) = record[1 + c];
        v(j, c) = record[4 + c];
        omega(j, c) = record[7 + c];
        f(j, c) = record[18 + c];
        tau(j, c) = record[21 + c];
      }
      for (unsigned int c = 0; c < 4; ++c)
        q(j, c) = record[10 + c];
      r(j) = record[14];
      m(j) = record[15];
      inertia(j) = record[16];
      hops(j) = record[17];
      frozen(j) = record[24];
      target_rank(j) = libMesh::DofObject::invalid_processor_id;
    }
  };

  /// Copy all particles of a host mirror into slots [begin, begin + mirror size) of this cloud
  void copyFrom(const HostMirror & host, const std::size_t begin) const
  {
    const auto range = std::make_pair(begin, begin + host.gid.extent(0));
    ::Kokkos::deep_copy(::Kokkos::subview(gid, range), host.gid);
    ::Kokkos::deep_copy(::Kokkos::subview(x, range, ::Kokkos::ALL), host.x);
    ::Kokkos::deep_copy(::Kokkos::subview(v, range, ::Kokkos::ALL), host.v);
    ::Kokkos::deep_copy(::Kokkos::subview(omega, range, ::Kokkos::ALL), host.omega);
    ::Kokkos::deep_copy(::Kokkos::subview(q, range, ::Kokkos::ALL), host.q);
    ::Kokkos::deep_copy(::Kokkos::subview(r, range), host.r);
    ::Kokkos::deep_copy(::Kokkos::subview(m, range), host.m);
    ::Kokkos::deep_copy(::Kokkos::subview(inertia, range), host.inertia);
    ::Kokkos::deep_copy(::Kokkos::subview(f, range, ::Kokkos::ALL), host.f);
    ::Kokkos::deep_copy(::Kokkos::subview(tau, range, ::Kokkos::ALL), host.tau);
    ::Kokkos::deep_copy(::Kokkos::subview(elem, range), host.elem);
    ::Kokkos::deep_copy(::Kokkos::subview(target_rank, range), host.target_rank);
    ::Kokkos::deep_copy(::Kokkos::subview(hops, range), host.hops);
    ::Kokkos::deep_copy(::Kokkos::subview(frozen, range), host.frozen);
  }

  /// Copy the first n particles of another cloud into the first n slots of this one
  void copyFrom(const ParticleCloud & other, const std::size_t n) const
  {
    const auto range = std::make_pair(std::size_t(0), n);
    ::Kokkos::deep_copy(::Kokkos::subview(gid, range), ::Kokkos::subview(other.gid, range));
    ::Kokkos::deep_copy(::Kokkos::subview(x, range, ::Kokkos::ALL),
                        ::Kokkos::subview(other.x, range, ::Kokkos::ALL));
    ::Kokkos::deep_copy(::Kokkos::subview(v, range, ::Kokkos::ALL),
                        ::Kokkos::subview(other.v, range, ::Kokkos::ALL));
    ::Kokkos::deep_copy(::Kokkos::subview(omega, range, ::Kokkos::ALL),
                        ::Kokkos::subview(other.omega, range, ::Kokkos::ALL));
    ::Kokkos::deep_copy(::Kokkos::subview(q, range, ::Kokkos::ALL),
                        ::Kokkos::subview(other.q, range, ::Kokkos::ALL));
    ::Kokkos::deep_copy(::Kokkos::subview(r, range), ::Kokkos::subview(other.r, range));
    ::Kokkos::deep_copy(::Kokkos::subview(m, range), ::Kokkos::subview(other.m, range));
    ::Kokkos::deep_copy(::Kokkos::subview(inertia, range),
                        ::Kokkos::subview(other.inertia, range));
    ::Kokkos::deep_copy(::Kokkos::subview(f, range, ::Kokkos::ALL),
                        ::Kokkos::subview(other.f, range, ::Kokkos::ALL));
    ::Kokkos::deep_copy(::Kokkos::subview(tau, range, ::Kokkos::ALL),
                        ::Kokkos::subview(other.tau, range, ::Kokkos::ALL));
    ::Kokkos::deep_copy(::Kokkos::subview(elem, range), ::Kokkos::subview(other.elem, range));
    ::Kokkos::deep_copy(::Kokkos::subview(target_rank, range),
                        ::Kokkos::subview(other.target_rank, range));
    ::Kokkos::deep_copy(::Kokkos::subview(hops, range), ::Kokkos::subview(other.hops, range));
    ::Kokkos::deep_copy(::Kokkos::subview(frozen, range), ::Kokkos::subview(other.frozen, range));
  }
};

} // namespace DEM
