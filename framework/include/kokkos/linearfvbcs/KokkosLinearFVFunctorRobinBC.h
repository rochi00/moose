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

/**
 * A Robin boundary condition alpha grad(u).n + beta u = gamma for Kokkos linear finite volume
 * kernels, with the three coefficients given as functors. It is the Kokkos analog of
 * LinearFVAdvectionDiffusionFunctorRobinBC, with the same expressions: the boundary value
 * eliminates the normal gradient with the normal distance to the face and extrapolates the cell
 * value with the cell gradient along the tangential part of the cell-to-face vector. As on the
 * host, the projections use the normal of the FaceInfo.
 */
class KokkosLinearFVFunctorRobinBC : public Moose::Kokkos::LinearFVBoundaryCondition
{
public:
  static InputParameters validParams();

  KokkosLinearFVFunctorRobinBC(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryValue(const FaceDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryNormalGradient(const FaceDatum & datum) const;

private:
  /// The coefficients of the normal gradient, of the value and the right-hand side
  ///@{
  Moose::Kokkos::LinearFVProperty _alpha;
  Moose::Kokkos::LinearFVProperty _beta;
  Moose::Kokkos::LinearFVProperty _gamma;
  ///@}
  /// The cell gradients of the variable
  Moose::Kokkos::LinearFVCellGradient _gradient;
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVFunctorRobinBC::BoundaryRelation
KokkosLinearFVFunctorRobinBC::computeBoundaryValue(const FaceDatum & datum) const
{
  // As LinearFVAdvectionDiffusionFunctorRobinBCBase::computeBoundaryValueMatrixContribution() and
  // computeBoundaryValueRHSContribution()
  const auto face = datum.face();
  const Real alpha = _alpha.face(face);
  const Real beta = _beta.face(face);
  const Real gamma = _gamma.face(face);
  const auto nhat = datum.faceInfoNormal();
  const auto d_cf = datum.cellToFaceVector();
  const auto grad_phi = _gradient(kokkosSystems(), datum.elem());
  const auto projection = d_cf * nhat;
  const auto vc = d_cf - (projection * nhat);
  return {alpha / (alpha + (beta * d_cf * nhat)),
          (gamma * projection / (alpha + (beta * projection))) +
              (alpha * grad_phi * vc / (alpha + (beta * projection)))};
}

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVFunctorRobinBC::BoundaryRelation
KokkosLinearFVFunctorRobinBC::computeBoundaryNormalGradient(const FaceDatum & datum) const
{
  // As LinearFVAdvectionDiffusionFunctorRobinBCBase::computeBoundaryGradientMatrixContribution()
  // and computeBoundaryGradientRHSContribution(); the matrix entry is the negated coefficient
  const auto face = datum.face();
  const Real alpha = _alpha.face(face);
  const Real beta = _beta.face(face);
  const Real gamma = _gamma.face(face);
  const auto nhat = datum.faceInfoNormal();
  const auto d_cf = datum.cellToFaceVector();
  const auto grad_phi = _gradient(kokkosSystems(), datum.elem());
  const auto projection = d_cf * nhat;
  const auto vc = d_cf - (projection * nhat);
  return {-(beta / (alpha + (beta * d_cf * nhat))),
          (gamma / alpha) + (-beta * gamma * projection / alpha / (alpha + (beta * projection))) +
              (-beta * grad_phi * vc / (alpha + (beta * projection)))};
}
