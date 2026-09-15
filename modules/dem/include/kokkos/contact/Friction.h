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
 * Coulomb friction on the tangential spring-dashpot of a contact model (plan layer L4). The
 * tangential force is -k_t delta_t - gamma_t v_t, limited to mu times the (repulsive part of
 * the) normal force; when it is, the spring is reset so that the spring and dashpot together
 * give exactly the limit, as in LAMMPS's granular pair styles.
 */
struct Friction
{
  /// Coefficient of friction; zero for a frictionless contact
  Real mu = 0;

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
    Moose::Kokkos::Real3 f_t = -k_t * delta_t - gamma_t * contact.v_t;
    const Real f_t_max = mu * (f_n > 0 ? f_n : 0);
    const Real f_t_norm = f_t.norm();
    if (f_t_norm > f_t_max)
    {
      f_t *= f_t_norm > 0 ? f_t_max / f_t_norm : 0;
      if (update)
        delta_t = (-1.0 / k_t) * (f_t + gamma_t * contact.v_t);
    }
    return f_t;
  }
};

} // namespace DEM
