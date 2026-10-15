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
 * Kokkos linear finite volume flux kernel implementing first-order upwind advection with a
 * constant velocity. It is the Kokkos analog of LinearFVAdvection for the upwind interpolation
 * case.
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
  /// Normal advective face flux
  KOKKOS_FUNCTION Real normalFaceFlux(const FaceDatum & datum) const;

  /// Constant advecting velocity vector
  const Moose::Kokkos::Real3 _velocity;
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVAdvection::InternalFaceContribution
KokkosLinearFVAdvection::computeInternalFaceContribution(const FaceDatum & datum) const
{
  const auto face_flux = normalFaceFlux(datum) * datum.faceArea();
  // Upwind: the outflow of the elem cell depends on the elem value, the inflow on the neighbor's
  return {face_flux > 0 ? face_flux : 0, face_flux < 0 ? face_flux : 0, 0, 0};
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
