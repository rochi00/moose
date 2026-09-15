#pragma once

#include "KokkosTypes.h"

namespace DEM
{

/**
 * Linear spring-dashpot normal contact model (plan layer L4), a compile-time policy so the force
 * kernel inlines it (decision D5). For two spheres of effective mass m_eff it is a damped linear
 * oscillator, so the coefficient of restitution is exp(-zeta pi / sqrt(1 - zeta^2)) with
 * zeta = damping / (2 sqrt(stiffness m_eff)) and the contact lasts pi / (omega_0 sqrt(1 - zeta^2)).
 */
struct LinearSpringDashpot
{
  /// Normal spring stiffness
  Real stiffness = 0;
  /// Normal dashpot coefficient
  Real damping = 0;

  /**
   * Normal force on particle i from particle j
   * @param overlap r_i + r_j minus the center distance, positive in contact
   * @param normal_velocity Relative velocity (v_i - v_j) projected on the unit normal from j to i;
   *        negative while approaching
   * @returns The signed force magnitude along the normal from j to i
   */
  KOKKOS_INLINE_FUNCTION Real normalForce(const Real overlap, const Real normal_velocity) const
  {
    return stiffness * overlap - damping * normal_velocity;
  }
};

} // namespace DEM
