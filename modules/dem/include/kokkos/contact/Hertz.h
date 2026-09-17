#pragma once

#include "KokkosTypes.h"

namespace DEM
{

/**
 * Hertz normal contact model with the Mindlin no-slip tangential stiffness (plan layer L4), a
 * compile-time policy like LinearSpringDashpot. The elastic force is
 *   F = (4/3) E* sqrt(r_eff) overlap^(3/2)
 * with E* the effective modulus E / (2 (1 - nu^2)) of two bodies of the same material (walls
 * included). The normal damping is one of
 *   Tsuji et al. (1992):  -2 sqrt(5/6) beta sqrt(S_n m_eff) v_n,   S_n = 2 E* sqrt(r_eff overlap),
 *     with beta = -ln e / sqrt(ln^2 e + pi^2) for a restitution coefficient e, which the contact
 *     then reproduces nearly independently of the impact velocity;
 *   Kuwabara and Kono (1987), Brilliantov et al. (1996):  A d(F_el)/d(overlap) d(overlap)/dt
 *     = -(3/2) A (4/3) E* sqrt(r_eff overlap) v_n, the viscoelastic damping with a dissipative
 *     constant A (a time), whose restitution coefficient falls with the impact velocity as
 *     1 - e ~ v^(1/5).
 */
struct Hertz
{
  /// Effective modulus E / (2 (1 - nu^2)); zero disables contact
  Real effective_modulus = 0;
  /// Damping ratio -ln e / sqrt(ln^2 e + pi^2), zero for an elastic contact
  Real beta = 0;
  /// Effective shear modulus G / (2 (2 - nu)) of two bodies of the same material, with
  /// G = E / (2 (1 + nu)); zero for a frictionless contact
  Real effective_shear_modulus = 0;
  /// Dissipative constant A of the Kuwabara-Kono damping, used instead of Tsuji's when positive
  Real dissipation_time = 0;
  /// Whether the damping may only reduce the repulsion, never make the normal force attractive
  /// (LAMMPS's limit_damping)
  bool limit_damping = false;

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
    const Real f_n =
        dissipation_time > 0
            ? 2.0 / 3.0 * s_n * overlap - dissipation_time * s_n * normal_velocity
            : 2.0 / 3.0 * s_n * overlap -
                  1.8257418583505538 * beta * ::Kokkos::sqrt(s_n * m_eff) * normal_velocity;
    return limit_damping && f_n < 0 ? 0 : f_n;
  }

  ///@{
  /// Mindlin no-slip tangential stiffness S_t = 8 G* sqrt(r_eff overlap) and the Tsuji damping
  /// -2 sqrt(5/6) beta sqrt(S_t m_eff) built on it
  KOKKOS_INLINE_FUNCTION Real tangentialStiffness(const Real overlap, const Real r_eff) const
  {
    return 8 * effective_shear_modulus * ::Kokkos::sqrt(r_eff * overlap);
  }
  KOKKOS_INLINE_FUNCTION Real tangentialDamping(const Real overlap,
                                                const Real r_eff,
                                                const Real m_eff) const
  {
    return 1.8257418583505538 * beta * ::Kokkos::sqrt(tangentialStiffness(overlap, r_eff) * m_eff);
  }
  ///@}

  /// Elastic energy at an overlap, the integral of the elastic force: (8/15) E* sqrt(r_eff) d^(5/2)
  KOKKOS_INLINE_FUNCTION Real energy(const Real overlap, const Real r_eff) const
  {
    return 8.0 / 15.0 * effective_modulus * ::Kokkos::sqrt(r_eff) * overlap * overlap *
           ::Kokkos::sqrt(overlap);
  }
};

} // namespace DEM
