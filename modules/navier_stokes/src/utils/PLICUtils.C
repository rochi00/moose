//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PLICUtils.h"
#include <cmath>
#include <algorithm>

namespace NS
{
namespace PLIC
{

Real
computePlaneOffset2DRect(const VectorValue<Real> & n_hat, Real alpha, Real dx, Real dy)
{
  // Clamp alpha to valid range
  alpha = std::max(0.0, std::min(1.0, alpha));

  if (alpha <= 0.0)
    return 0.0;
  if (alpha >= 1.0)
    return std::abs(n_hat(0)) * dx + std::abs(n_hat(1)) * dy;

  // Work with absolute values of normal components scaled by cell dimensions
  const Real m1_raw = std::abs(n_hat(0)) * dx;
  const Real m2_raw = std::abs(n_hat(1)) * dy;

  // Ensure m1 <= m2
  const Real m1 = std::min(m1_raw, m2_raw);
  const Real m2 = std::max(m1_raw, m2_raw);

  // Handle degenerate normal (nearly aligned with one axis)
  if (m1 < 1.0e-14)
    return alpha * m2;

  const Real m12 = m1 * m2;

  // Three regimes based on alpha thresholds
  // Threshold between triangle and trapezoid: alpha_1 = m1^2 / (2 * m1 * m2) = m1 / (2 * m2)
  const Real alpha_1 = m1 / (2.0 * m2);
  // Threshold between trapezoid and complement: alpha_2 = 1 - m1 / (2 * m2)
  const Real alpha_2 = 1.0 - alpha_1;

  Real d;
  if (alpha <= alpha_1)
    // Triangle regime
    d = std::sqrt(2.0 * m12 * alpha);
  else if (alpha <= alpha_2)
    // Trapezoid regime
    d = alpha * m2 + 0.5 * m1;
  else
    // Complement triangle regime
    d = (m1 + m2) - std::sqrt(2.0 * m12 * (1.0 - alpha));

  return d;
}

Real
truncatedVolumeFraction2DRect(const VectorValue<Real> & n_hat,
                              Real d_local,
                              Real dx,
                              Real dy)
{
  const Real m1_raw = std::abs(n_hat(0)) * dx;
  const Real m2_raw = std::abs(n_hat(1)) * dy;
  const Real m1 = std::min(m1_raw, m2_raw);
  const Real m2 = std::max(m1_raw, m2_raw);
  const Real m_sum = m1 + m2;
  const Real m12 = m1 * m2;

  // Clamp d to valid range
  d_local = std::max(0.0, std::min(d_local, m_sum));

  if (m12 < 1.0e-28)
    return (m2 > 0.0) ? d_local / m2 : 0.0;

  Real alpha;
  if (d_local <= m1)
    alpha = d_local * d_local / (2.0 * m12);
  else if (d_local <= m2)
    alpha = (d_local - 0.5 * m1) / m2;
  else
    alpha = 1.0 - (m_sum - d_local) * (m_sum - d_local) / (2.0 * m12);

  return std::max(0.0, std::min(1.0, alpha));
}

Real
signedDistanceToPlane(const Point & p, const VectorValue<Real> & n_hat, Real d)
{
  return n_hat * p - d;
}

Real
localToGlobalOffset(const VectorValue<Real> & n_hat,
                    Real d_local,
                    const Point & x0,
                    Real dx,
                    Real dy)
{
  // The local frame uses |n_x|, |n_y| with origin at the cell corner
  // where n_hat . x is minimized. For each component:
  //   if n_i >= 0: the minimum is at x0_i, contribute n_i * x0_i
  //   if n_i <  0: the minimum is at x0_i + L_i, contribute n_i * (x0_i + L_i)
  Real d_global = d_local;
  d_global += (n_hat(0) >= 0.0) ? n_hat(0) * x0(0) : n_hat(0) * (x0(0) + dx);
  d_global += (n_hat(1) >= 0.0) ? n_hat(1) * x0(1) : n_hat(1) * (x0(1) + dy);
  return d_global;
}

Real
globalToLocalOffset(const VectorValue<Real> & n_hat,
                    Real d_global,
                    const Point & x0,
                    Real dx,
                    Real dy)
{
  // Reverse of localToGlobalOffset: subtract the origin contribution
  Real d_local = d_global;
  d_local -= (n_hat(0) >= 0.0) ? n_hat(0) * x0(0) : n_hat(0) * (x0(0) + dx);
  d_local -= (n_hat(1) >= 0.0) ? n_hat(1) * x0(1) : n_hat(1) * (x0(1) + dy);
  return d_local;
}

Real
volumeFractionInRect(const VectorValue<Real> & n_hat,
                     Real d_global,
                     const Point & rect_origin,
                     Real rect_dx,
                     Real rect_dy)
{
  if (rect_dx < 1.0e-30 || rect_dy < 1.0e-30)
    return 0.0;

  const Real d_local = globalToLocalOffset(n_hat, d_global, rect_origin, rect_dx, rect_dy);
  return truncatedVolumeFraction2DRect(n_hat, d_local, rect_dx, rect_dy);
}

} // namespace PLIC
} // namespace NS
