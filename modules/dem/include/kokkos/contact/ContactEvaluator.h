#pragma once

#include "ParticleCloud.h"
#include "AnalyticWalls.h"
#include "SidesetWalls.h"
#include "PairState.h"
#include "Friction.h"

namespace DEM
{

/**
 * Device-side evaluation of the contacts of a cloud with one contact model: the sphere-sphere
 * and sphere-wall narrow phases, and the normal and tangential forces of a contact from its
 * history in the pair-state map. Shared by the kernel that advances the histories once per
 * contact and substep and the kernel that sums the forces on every particle from both sides.
 */
template <typename Model>
struct ContactEvaluator
{
  ParticleCloud cloud;
  AnalyticWalls walls;
  SidesetWalls sidesets;
  Model model;
  /// Friction of particle-particle and of particle-wall contacts
  Friction friction;
  Friction wall_friction;
  PairStateMap states;
  /// Substep, by which the tangential springs are stretched
  Real dt;

  /// The contact of particle i with particle j as seen from i, if they overlap
  KOKKOS_INLINE_FUNCTION bool sphere(const std::size_t i, const std::size_t j, Contact & c) const
  {
    const Moose::Kokkos::Real3 d(
        cloud.x(i, 0) - cloud.x(j, 0), cloud.x(i, 1) - cloud.x(j, 1), cloud.x(i, 2) - cloud.x(j, 2));
    const Real distance = d.norm();
    const Real overlap = cloud.r(i) + cloud.r(j) - distance;
    if (overlap <= 0)
      return false;
    c = Contact((1.0 / distance) * d,
                overlap,
                Moose::Kokkos::Real3(cloud.v(i, 0), cloud.v(i, 1), cloud.v(i, 2)),
                Moose::Kokkos::Real3(cloud.omega(i, 0), cloud.omega(i, 1), cloud.omega(i, 2)),
                cloud.r(i),
                cloud.m(i),
                Moose::Kokkos::Real3(cloud.v(j, 0), cloud.v(j, 1), cloud.v(j, 2)),
                Moose::Kokkos::Real3(cloud.omega(j, 0), cloud.omega(j, 1), cloud.omega(j, 2)),
                cloud.r(j),
                cloud.m(j));
    return true;
  }

  /// The contacts of particle i with the sideset walls within a reach of its center, reduced
  /// to the physical ones or not (the faces to list histories for)
  KOKKOS_INLINE_FUNCTION unsigned int sidesetContacts(const std::size_t i,
                                                      const Real reach,
                                                      WallContact * const out,
                                                      const bool reduce = true) const
  {
    return sidesets.contacts(cloud.elem(i),
                             Moose::Kokkos::Real3(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2)),
                             cloud.r(i),
                             reach,
                             out,
                             reduce);
  }

  /// The contact of particle i with a sideset wall face found by sidesetContacts()
  KOKKOS_INLINE_FUNCTION Contact sideset(const std::size_t i, const WallContact & wc) const
  {
    return Contact(wc.normal,
                   wc.overlap,
                   Moose::Kokkos::Real3(cloud.v(i, 0), cloud.v(i, 1), cloud.v(i, 2)),
                   Moose::Kokkos::Real3(cloud.omega(i, 0), cloud.omega(i, 1), cloud.omega(i, 2)),
                   cloud.r(i),
                   cloud.m(i),
                   wc.velocity,
                   Moose::Kokkos::Real3(0),
                   0,
                   0);
  }

  /// History key of particle i's contact with a sideset wall face, numbered after the analytic
  /// walls by global face ID
  KOKKOS_INLINE_FUNCTION PairKey sidesetKey(const std::size_t i, const int64_t id) const
  {
    return PairKey::wall(cloud.gid(i), walls.n + static_cast<std::size_t>(id));
  }
  KOKKOS_INLINE_FUNCTION PairKey sidesetKey(const std::size_t i, const WallContact & wc) const
  {
    return sidesetKey(i, wc.id);
  }

  /**
   * Carry the history of a contact that slid across an edge onto its new face: for the touching
   * contact k with an unstretched history, a neighbor face on the same surface that is not
   * among the touching contacts and has a stretched history hands it over and is reset
   * @param wc The contacts of particle i, the first n_touching of which overlap
   */
  KOKKOS_INLINE_FUNCTION void transferHistory(const std::size_t i,
                                              const WallContact * const wc,
                                              const unsigned int n_touching,
                                              const unsigned int k) const
  {
    const auto to = states.find(sidesetKey(i, wc[k]));
    if (!states.valid_at(to) || states.value_at(to).stretched())
      return;
    for (unsigned int e = 0; e < 3; ++e)
    {
      const int g = sidesets.neighbors(wc[k].face, e);
      if (g < 0 || !sidesets.same_surface(wc[k].face, e))
        continue;
      bool touching = false;
      for (unsigned int l = 0; l < n_touching && !touching; ++l)
        touching = wc[l].face == static_cast<std::size_t>(g);
      if (touching)
        continue;
      const auto from = states.find(sidesetKey(i, sidesets.ids(g)));
      if (!states.valid_at(from) || !states.value_at(from).stretched())
        continue;
      states.value_at(to) = states.value_at(from);
      states.value_at(from) = PairState();
      return;
    }
  }

  /// The contact of particle i with wall w, if they overlap
  KOKKOS_INLINE_FUNCTION bool wall(const std::size_t i, const std::size_t w, Contact & c) const
  {
    const Moose::Kokkos::Real3 xi(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2));
    const Real overlap = walls.overlap(w, xi, cloud.r(i));
    if (overlap <= 0)
      return false;
    c = Contact(walls.normal(w),
                overlap,
                Moose::Kokkos::Real3(cloud.v(i, 0), cloud.v(i, 1), cloud.v(i, 2)),
                Moose::Kokkos::Real3(cloud.omega(i, 0), cloud.omega(i, 1), cloud.omega(i, 2)),
                cloud.r(i),
                cloud.m(i),
                walls.velocity(w),
                Moose::Kokkos::Real3(0),
                0,
                0);
    return true;
  }

  /**
   * Normal and tangential forces and rolling resistance torque on the particle a contact is
   * seen from
   * @param c The contact
   * @param key Its history key
   * @param sign +1 when the particle is the key's lower one (or the contact is with a wall), -1
   *        otherwise: the tangential history is stored in the lower particle's orientation
   * @param update Whether to advance the history by the substep or only read it
   */
  KOKKOS_INLINE_FUNCTION void forces(const Contact & c,
                                     const PairKey & key,
                                     const Real sign,
                                     const bool update,
                                     Real & f_n,
                                     Moose::Kokkos::Real3 & f_t,
                                     Moose::Kokkos::Real3 & tau_r) const
  {
    f_n = model.normalForce(c.overlap, c.v_n, c.r_eff, c.m_eff);
    // A contact without a history entry (which the rebuild gives every listed pair) is treated
    // as fresh
    const auto entry = states.find(key);
    const bool found = states.valid_at(entry);
    Moose::Kokkos::Real3 delta_t(0), delta_r(0);
    if (found)
      for (unsigned int k = 0; k < 3; ++k)
      {
        delta_t(k) = sign * states.value_at(entry).delta_t[k];
        delta_r(k) = states.value_at(entry).delta_r[k];
      }
    // Wall keys carry the negated wall index as their upper entry
    const Friction & fr = key.upper < 0 ? wall_friction : friction;
    f_t = fr.tangentialForce(model, c, f_n, dt, delta_t, update);
    tau_r = fr.rollingTorque(c, f_n, dt, delta_r, update);
    if (update && found)
      for (unsigned int k = 0; k < 3; ++k)
      {
        states.value_at(entry).delta_t[k] = sign * delta_t(k);
        states.value_at(entry).delta_r[k] = delta_r(k);
      }
  }

  /// Reset the history of a separated contact
  KOKKOS_INLINE_FUNCTION void reset(const PairKey & key) const
  {
    const auto entry = states.find(key);
    if (states.valid_at(entry))
      states.value_at(entry) = PairState();
  }
};

} // namespace DEM
