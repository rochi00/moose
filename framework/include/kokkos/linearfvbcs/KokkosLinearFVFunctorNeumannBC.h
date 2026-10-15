//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosLinearFVBoundaryCondition.h"
#include "KokkosParsedFunction.h"

/**
 * Prescribes the outward diffusive flux on a boundary for Kokkos linear finite volume kernels. It is
 * the Kokkos analog of LinearFVAdvectionDiffusionFunctorNeumannBC: the boundary value extrapolates
 * the cell value with the flux divided by the diffusion coefficient along the normal, plus the cell
 * gradient along the tangential part of the cell-to-face vector.
 */
class KokkosLinearFVFunctorNeumannBC : public Moose::Kokkos::LinearFVBoundaryCondition
{
public:
  static InputParameters validParams();

  KokkosLinearFVFunctorNeumannBC(const InputParameters & parameters);

  virtual bool includesMaterialPropertyMultiplier() const override { return true; }

  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryValue(const FaceDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryNormalGradient(const FaceDatum & datum) const;

private:
  /// The functor providing the outward diffusive flux on the boundary
  const Moose::Kokkos::ReferenceWrapper<const KokkosParsedFunction> _flux_functor;
  /// The diffusion coefficient
  const Moose::Kokkos::ReferenceWrapper<const KokkosParsedFunction> _diffusion_coeff;
  /// The cell gradients of the variable
  Moose::Kokkos::LinearFVCellGradient _gradient;
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVFunctorNeumannBC::BoundaryRelation
KokkosLinearFVFunctorNeumannBC::computeBoundaryValue(const FaceDatum & datum) const
{
  const auto face_centroid = datum.faceCentroid();
  const auto normal = datum.faceNormal();
  const auto d_cf = datum.cellToFaceVector();
  // The part of the cell-to-face vector orthogonal to the normal, along which the cell gradient
  // extrapolates on non-orthogonal meshes
  const auto correction_vector = d_cf - (d_cf * normal) * normal;
  return {1,
          _flux_functor->value(_t, face_centroid) / _diffusion_coeff->value(_t, face_centroid) *
                  datum.cellToFaceDistance() +
              _gradient(kokkosSystems(), datum.elem()) * correction_vector};
}

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVFunctorNeumannBC::BoundaryRelation
KokkosLinearFVFunctorNeumannBC::computeBoundaryNormalGradient(const FaceDatum & datum) const
{
  // A diffusive flux, see includesMaterialPropertyMultiplier()
  return {0, _flux_functor->value(_t, datum.faceCentroid())};
}
