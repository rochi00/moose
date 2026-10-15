//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosLinearFVFluxKernel.h"

/**
 * Kokkos linear finite volume flux kernel implementing the advection term with a constant
 * velocity. It is the Kokkos analog of LinearFVAdvection: the advected face value on internal
 * faces comes from the FVInterpolationMethod named by advected_interp_method_name, which must be
 * one of FVAdvectedUpwind, FVGeometricAverage, FVAdvectedMinmodWeightBased,
 * FVAdvectedVanLeerWeightBased or FVAdvectedMUSCLDeferredCorrection; its parameters are copied to
 * the device.
 */
class KokkosLinearFVAdvection : public Moose::Kokkos::LinearFVFluxKernel
{
public:
  static InputParameters validParams();

  KokkosLinearFVAdvection(const InputParameters & parameters);

  virtual bool needsBoundaryValueData() const override { return true; }

  template <typename Derived>
  KOKKOS_FUNCTION InternalFaceContribution
  computeInternalFaceContribution(const FaceDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryFaceContribution computeBoundaryFaceContribution(
      const FaceDatum & datum, const int bc_index, const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const;

private:
  /// The supported advected interpolation methods
  enum class Interpolation
  {
    UPWIND,
    AVERAGE,
    MINMOD,
    VAN_LEER,
    MUSCL
  };

  /// Normal advective face flux
  KOKKOS_FUNCTION Real normalFaceFlux(const FaceDatum & datum) const;

  /// Constant advecting velocity vector
  const Moose::Kokkos::Real3 _velocity;
  /// The advected interpolation method
  Interpolation _interpolation;
  /// FVAdvectedMinmodWeightBased and FVAdvectedVanLeerWeightBased: whether the weights are no more
  /// downwind-biased than linear interpolation
  bool _limit_to_linear = true;
  /// FVAdvectedMinmodWeightBased and FVAdvectedVanLeerWeightBased: the high-order blending strength
  Real _blending_factor = 1;
  /// FVAdvectedMUSCLDeferredCorrection: the deferred correction strength
  Real _deferred_correction_factor = 1;
  /// The cell values of the variable, for the interpolation methods that need them
  const Moose::Kokkos::LinearFVCellValues _u;
  /// The cell gradients of the variable from the interpolation method's gradient method, for the
  /// interpolation methods that need them
  Moose::Kokkos::LinearFVCellGradient _gradient;
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVAdvection::InternalFaceContribution
KokkosLinearFVAdvection::computeInternalFaceContribution(const FaceDatum & datum) const
{
  const auto mass_flux = normalFaceFlux(datum);
  const auto area = datum.faceArea();

  // Weights of the elem and neighbor values in the face value and the explicit part of the face
  // value, which the right-hand side carries; as FVAdvectedInterpolationMethod::advectedInterpolate()
  Real w_elem, w_neighbor, rhs_face_value = 0;
  if (_interpolation == Interpolation::UPWIND)
  {
    const Real neighbor_weight = mass_flux < 0.0;
    w_elem = 1.0 - neighbor_weight;
    w_neighbor = neighbor_weight;
  }
  else if (_interpolation == Interpolation::AVERAGE)
  {
    const auto gc = datum.faceGC();
    w_elem = gc;
    w_neighbor = 1.0 - gc;
  }
  else
  {
    const auto & systems = kokkosSystems();
    const auto elem = datum.elem();
    const auto neighbor = datum.neighbor();
    const bool upwind_is_elem = mass_flux >= 0.0;
    const auto elem_value = _u(systems, elem);
    const auto neighbor_value = _u(systems, neighbor);
    const auto phi_upwind = upwind_is_elem ? elem_value : neighbor_value;
    const auto grad_upwind = _gradient(systems, upwind_is_elem ? elem : neighbor);

    if (_interpolation == Interpolation::MUSCL)
    {
      // Matrix: upwind; right-hand side: the difference to the face value reconstructed from the
      // upwind cell with its gradient
      const auto phi_high =
          phi_upwind + grad_upwind * (datum.faceCentroid() - (upwind_is_elem
                                                                  ? datum.elemCentroid()
                                                                  : datum.neighborCentroid()));
      w_elem = upwind_is_elem ? 1.0 : 0.0;
      w_neighbor = upwind_is_elem ? 0.0 : 1.0;
      const auto phi_matrix = w_elem * elem_value + w_neighbor * neighbor_value;
      rhs_face_value = _deferred_correction_factor * (phi_matrix - phi_high);
    }
    else
    {
      // Limited blending of the upwind and downwind values, as Moose::FV::rF() and the weight
      // based Minmod and Van Leer methods
      const auto phi_downwind = upwind_is_elem ? neighbor_value : elem_value;
      const auto d_cn = datum.neighborCentroid() - datum.elemCentroid();
      const auto grad_dot = grad_upwind * (upwind_is_elem ? d_cn : -d_cn);
      const auto denom = phi_downwind - phi_upwind;
      // The tolerance of Moose::FV::rF()
      constexpr Real eps = 1e-10;
      const Real denom_sign =
          ::Kokkos::copysign(1.0, denom + std::numeric_limits<Real>::min() * grad_dot);
      const auto r_f = 2.0 * grad_dot / (denom + denom_sign * eps) - 1.0;
      const auto abs_r_f = r_f < 0 ? -r_f : r_f;
      const Real beta = _interpolation == Interpolation::MINMOD
                            ? (r_f < 1.0 ? (r_f > 0.0 ? r_f : 0.0) : 1.0)
                            : (r_f + abs_r_f) / (1.0 + abs_r_f);

      // Greenshields blending: phi_f = (1 - g) phi_upwind + g phi_downwind, g = beta (1 - w_f),
      // with w_f the geometric weight of the upwind cell
      const auto gc = datum.faceGC();
      const Real w_f = upwind_is_elem ? gc : (1.0 - gc);
      const Real g_unclamped = _blending_factor * beta * (1.0 - w_f);
      const Real g_lower = g_unclamped > 0.0 ? g_unclamped : 0.0;
      const Real g_clamped = g_lower < 1.0 - w_f ? g_lower : 1.0 - w_f;
      const Real g = _limit_to_linear ? g_clamped : g_unclamped;
      w_elem = upwind_is_elem ? 1.0 - g : g;
      w_neighbor = upwind_is_elem ? g : 1.0 - g;
    }
  }

  const auto rhs = rhs_face_value * mass_flux * area;
  return {w_elem * mass_flux * area, w_neighbor * mass_flux * area, rhs, -rhs};
}

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVAdvection::BoundaryFaceContribution
KokkosLinearFVAdvection::computeBoundaryFaceContribution(
    const FaceDatum & datum, const int bc_index, const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const
{
  const auto face_flux = normalFaceFlux(datum) * datum.faceArea();
  return {face_flux * boundaryValueCoefficient(bc_index, bc_face_index),
          -face_flux * boundaryValueSource(bc_index, bc_face_index)};
}

KOKKOS_FUNCTION inline Real
KokkosLinearFVAdvection::normalFaceFlux(const FaceDatum & datum) const
{
  return _velocity * datum.faceNormal();
}
