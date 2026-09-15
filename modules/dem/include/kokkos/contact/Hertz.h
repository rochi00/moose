#pragma once

#include "KokkosTypes.h"

namespace DEM
{

/**
 * Hertz normal contact model with the damping of Tsuji et al. (1992) (plan layer L4), a
 * compile-time policy like LinearSpringDashpot. The elastic force is
 *   F = (4/3) E* sqrt(r_eff) overlap^(3/2)
 * with E* the effective modulus E / (2 (1 - nu^2)) of two bodies of the same material (walls
 * included), and the damping force is
 *   -2 sqrt(5/6) beta sqrt(S_n m_eff) v_n,   S_n = 2 E* sqrt(r_eff overlap),
 * with beta = -ln e / sqrt(ln^2 e + pi^2) for a restitution coefficient e, which the contact
 * then reproduces nearly independently of the impact velocity.
 */
struct Hertz
{
  /// Effective modulus E / (2 (1 - nu^2)); zero disables contact
  Real effective_modulus = 0;
  /// Damping ratio -ln e / sqrt(ln^2 e + pi^2), zero for an elastic contact
  Real beta = 0;

  /// Whether the model applies any force
  bool enabled() const { return effective_modulus > 0; }

  /// Normal force on particle i from particle j; see LinearSpringDashpot::normalForce
  KOKKOS_INLINE_FUNCTION Real normalForce(const Real overlap,
                                          const Real normal_velocity,
                                          const Real r_eff,
                                          const Real m_eff) const
  {
    // S_n = 2 E* sqrt(r_eff overlap) is the tangent stiffness dF/d(overlap)
    const Real s_n = 2 * effective_modulus * ::Kokkos::sqrt(r_eff * overlap);
    // 2 sqrt(5/6) = 1.8257418583505538
    return 2.0 / 3.0 * s_n * overlap -
           1.8257418583505538 * beta * ::Kokkos::sqrt(s_n * m_eff) * normal_velocity;
  }

  /// Elastic energy at an overlap, the integral of the elastic force: (8/15) E* sqrt(r_eff) d^(5/2)
  KOKKOS_INLINE_FUNCTION Real energy(const Real overlap, const Real r_eff) const
  {
    return 8.0 / 15.0 * effective_modulus * ::Kokkos::sqrt(r_eff) * overlap * overlap *
           ::Kokkos::sqrt(overlap);
  }
};

} // namespace DEM
