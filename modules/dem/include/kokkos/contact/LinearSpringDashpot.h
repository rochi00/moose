#pragma once

#include "KokkosTypes.h"

namespace DEM
{

/**
 * Linear spring-dashpot normal contact model (plan layer L4), a compile-time policy so the force
 * kernel inlines it (decision D5). For two spheres of effective mass m_eff it is a damped linear
 * oscillator, so the coefficient of restitution is exp(-zeta pi / sqrt(1 - zeta^2)) with
 * zeta = damping / (2 sqrt(stiffness m_eff)) and the contact lasts pi / (omega_0 sqrt(1 - zeta^2)).
 *
 * Every contact model provides enabled(), normalForce(), energy(), tangentialStiffness(), and
 * tangentialDamping() with the signatures below; the effective radius and mass of the pair are
 * r_i r_j / (r_i + r_j) and m_i m_j / (m_i + m_j), or r_i and m_i against a wall. The tangential
 * spring acts on the accumulated tangential displacement of the contact (PairState) under the
 * Coulomb limit applied by Friction.
 */
struct LinearSpringDashpot
{
  /// Normal spring stiffness
  Real stiffness = 0;
  /// Normal dashpot coefficient
  Real damping = 0;
  /// Tangential spring stiffness; zero for a frictionless contact
  Real tangential_stiffness = 0;
  /// Tangential dashpot coefficient
  Real tangential_damping = 0;
  /// Whether the dashpot may only reduce the repulsion, never make the normal force attractive
  /// (LAMMPS's limit_damping)
  bool limit_damping = false;

  /// Whether the model applies any force
  bool enabled() const { return stiffness > 0; }

  /**
   * Normal force on particle i from particle j
   * @param overlap r_i + r_j minus the center distance, positive in contact
   * @param normal_velocity Relative velocity (v_i - v_j) projected on the unit normal from j to i;
   *        negative while approaching
   * @param r_eff Effective radius of the pair (unused)
   * @param m_eff Effective mass of the pair (unused)
   * @returns The signed force magnitude along the normal from j to i
   */
  KOKKOS_INLINE_FUNCTION Real normalForce(const Real overlap,
                                          const Real normal_velocity,
                                          const Real /*r_eff*/,
                                          const Real /*m_eff*/) const
  {
    const Real f_n = stiffness * overlap - damping * normal_velocity;
    return limit_damping && f_n < 0 ? 0 : f_n;
  }

  /// Elastic energy stored in the spring at an overlap
  KOKKOS_INLINE_FUNCTION Real energy(const Real overlap, const Real /*r_eff*/) const
  {
    return 0.5 * stiffness * overlap * overlap;
  }

  ///@{
  /// Tangential spring stiffness and dashpot coefficient of a contact
  KOKKOS_INLINE_FUNCTION Real tangentialStiffness(const Real /*overlap*/,
                                                  const Real /*r_eff*/) const
  {
    return tangential_stiffness;
  }
  KOKKOS_INLINE_FUNCTION Real tangentialDamping(const Real /*overlap*/,
                                                const Real /*r_eff*/,
                                                const Real /*m_eff*/) const
  {
    return tangential_damping;
  }
  ///@}
};

} // namespace DEM
