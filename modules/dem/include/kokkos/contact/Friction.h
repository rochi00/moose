#pragma once

#include "KokkosTypes.h"

namespace DEM
{

/**
 * Geometry and kinematics of one contact as seen from particle i: the unit normal from the
 * partner toward i, the overlap, the relative velocity of the contact point on i with respect
 * to the contact point on the partner, split into its normal and tangential parts, and the
 * effective radius and mass of the pair
 */
struct Contact
{
  Moose::Kokkos::Real3 normal;
  Real overlap = 0;
  Moose::Kokkos::Real3 v_rel;
  Real v_n = 0;
  Moose::Kokkos::Real3 v_t;
  /// Relative spin omega_i - omega_j, for rolling resistance
  Moose::Kokkos::Real3 omega_rel;
  Real r_eff = 0;
  Real m_eff = 0;

  KOKKOS_INLINE_FUNCTION Contact() = default;

  /**
   * Contact of sphere i with a partner: a sphere j, or a wall (rigid, at rest, of infinite
   * radius and mass) when r_j is zero
   */
  KOKKOS_INLINE_FUNCTION Contact(const Moose::Kokkos::Real3 & normal,
                                 const Real overlap,
                                 const Moose::Kokkos::Real3 & v_i,
                                 const Moose::Kokkos::Real3 & omega_i,
                                 const Real r_i,
                                 const Real m_i,
                                 const Moose::Kokkos::Real3 & v_j,
                                 const Moose::Kokkos::Real3 & omega_j,
                                 const Real r_j,
                                 const Real m_j)
    : normal(normal), overlap(overlap)
  {
    // The contact point is at -r_i n from x_i and at r_j n from x_j, so its velocity on each
    // body picks up the spin; the surfaces are taken at the undeformed radii
    v_rel = v_i - v_j - (r_i * omega_i + r_j * omega_j).cross_product(normal);
    v_n = v_rel.dot_product(normal);
    v_t = v_rel - v_n * normal;
    omega_rel = omega_i - omega_j;
    if (r_j > 0)
    {
      r_eff = r_i * r_j / (r_i + r_j);
      m_eff = m_i * m_j / (m_i + m_j);
    }
    else
    {
      r_eff = r_i;
      m_eff = m_i;
    }
  }
};

/**
 * Coulomb friction on the tangential spring-dashpot of a contact model, and elastic-plastic
 * spring-dashpot rolling resistance (plan layer L4). The tangential force is
 * -k_t delta_t - gamma_t v_t, limited to mu times the (repulsive part of the) normal force; when
 * it is, the spring is reset so that the spring and dashpot together give exactly the limit, as
 * in LAMMPS's granular pair styles. With partial slip, the spring instead follows the
 * Mindlin-Deresiewicz loading curve mu F_n (1 - (1 - |delta_t| / delta_max)^(3/2)) with
 * delta_max = 3 mu F_n / (2 k_t), whose initial slope is k_t and which reaches the limit at
 * delta_max (the closed form of Di Renzo and Di Maio 2004, elastic on unloading; with the Hertz
 * model's k_t = 8 G* a it is Mindlin's solution for a constant normal force). Rolling
 * resistance works the same way as the linear spring on the rolling
 * displacement, the integral of the rolling velocity r_eff (omega_i - omega_j) x n, with its own
 * stiffness, damping, and coefficient, and acts as the torque r_eff n x F_r on i and its
 * opposite on j (Luding 2008; LAMMPS's rolling sds).
 */
struct Friction
{
  /// Coefficient of friction; zero for a frictionless contact
  Real mu = 0;
  /// Whether the tangential spring follows the Mindlin-Deresiewicz partial-slip loading curve
  /// rather than the linear no-slip spring
  bool partial_slip = false;
  /// Coefficient of rolling friction, limiting the rolling torque to this times the normal force
  /// times the effective radius; zero for no rolling resistance
  Real mu_r = 0;
  /// Rolling spring stiffness and dashpot coefficient, on the rolling displacement
  Real k_r = 0;
  Real gamma_r = 0;

  /**
   * Tangential force on particle i of a contact
   * @param model The contact model, for the tangential stiffness and damping
   * @param contact The contact as seen from i
   * @param f_n The normal force on i (signed along the normal)
   * @param dt The substep, by which the spring is stretched when updating
   * @param delta_t The tangential spring displacement in the orientation of i; when updating,
   *        it is first projected on the current tangent plane and stretched by v_t dt, and left
   *        holding the value the force was computed from
   * @param update Whether to advance the spring or only evaluate the force from it
   */
  template <typename Model>
  KOKKOS_INLINE_FUNCTION Moose::Kokkos::Real3 tangentialForce(const Model & model,
                                                              const Contact & contact,
                                                              const Real f_n,
                                                              const Real dt,
                                                              Moose::Kokkos::Real3 & delta_t,
                                                              const bool update) const
  {
    const Real k_t = model.tangentialStiffness(contact.overlap, contact.r_eff);
    if (k_t <= 0 || mu <= 0)
      return Moose::Kokkos::Real3(0);
    const Real gamma_t = model.tangentialDamping(contact.overlap, contact.r_eff, contact.m_eff);
    if (update)
    {
      delta_t -= delta_t.dot_product(contact.normal) * contact.normal;
      delta_t += dt * contact.v_t;
    }
    const Real f_t_max = mu * (f_n > 0 ? f_n : 0);
    Moose::Kokkos::Real3 f_t;
    if (partial_slip)
    {
      // The spring saturates at delta_max: beyond it the contact slides and the spring is held
      // there, so unloading starts from the fully slipped state
      const Real delta_max = 1.5 * f_t_max / k_t;
      const Real delta_norm = delta_t.norm();
      if (delta_norm >= delta_max)
      {
        if (update && delta_norm > 0)
          delta_t *= delta_max / delta_norm;
        f_t = delta_norm > 0 ? (-f_t_max / delta_norm) * delta_t : Moose::Kokkos::Real3(0);
      }
      else
      {
        const Real remaining = 1 - delta_norm / delta_max;
        const Real spring = f_t_max * (1 - remaining * ::Kokkos::sqrt(remaining));
        f_t = delta_norm > 0 ? (-spring / delta_norm) * delta_t : Moose::Kokkos::Real3(0);
      }
      f_t -= gamma_t * contact.v_t;
      const Real f_t_norm = f_t.norm();
      if (f_t_norm > f_t_max)
        f_t *= f_t_max / f_t_norm;
      return f_t;
    }
    f_t = -k_t * delta_t - gamma_t * contact.v_t;
    const Real f_t_norm = f_t.norm();
    if (f_t_norm > f_t_max)
    {
      f_t *= f_t_norm > 0 ? f_t_max / f_t_norm : 0;
      if (update)
        delta_t = (-1.0 / k_t) * (f_t + gamma_t * contact.v_t);
    }
    return f_t;
  }

  /**
   * Rolling resistance torque on particle i of a contact
   * @param delta_r The rolling spring displacement, the same from either side; when updating, it
   *        is first projected on the tangent plane and stretched by the rolling velocity times
   *        dt, and left holding the value the torque was computed from
   */
  KOKKOS_INLINE_FUNCTION Moose::Kokkos::Real3 rollingTorque(const Contact & contact,
                                                            const Real f_n,
                                                            const Real dt,
                                                            Moose::Kokkos::Real3 & delta_r,
                                                            const bool update) const
  {
    if (k_r <= 0 || mu_r <= 0)
      return Moose::Kokkos::Real3(0);
    const Moose::Kokkos::Real3 u = contact.r_eff * contact.omega_rel.cross_product(contact.normal);
    if (update)
    {
      delta_r -= delta_r.dot_product(contact.normal) * contact.normal;
      delta_r += dt * u;
    }
    Moose::Kokkos::Real3 f_r = -k_r * delta_r - gamma_r * u;
    const Real f_r_max = mu_r * (f_n > 0 ? f_n : 0);
    const Real f_r_norm = f_r.norm();
    if (f_r_norm > f_r_max)
    {
      f_r *= f_r_norm > 0 ? f_r_max / f_r_norm : 0;
      if (update)
        delta_r = (-1.0 / k_r) * (f_r + gamma_r * u);
    }
    return contact.r_eff * contact.normal.cross_product(f_r);
  }
};

} // namespace DEM
