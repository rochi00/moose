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
 * A symmetry boundary for a scalar quantity of Kokkos linear finite volume kernels: no normal
 * gradient, and the boundary value is the cell value, optionally with the cell gradient along the
 * tangential part of the cell-to-face vector. It is the Kokkos analog of
 * LinearFVAdvectionDiffusionScalarSymmetryBC.
 */
class KokkosLinearFVScalarSymmetryBC : public Moose::Kokkos::LinearFVBoundaryCondition
{
public:
  static InputParameters validParams();

  KokkosLinearFVScalarSymmetryBC(const InputParameters & parameters);

  virtual bool includesMaterialPropertyMultiplier() const override { return true; }

  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryValue(const FaceDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryNormalGradient(const FaceDatum &) const
  {
    return {0, 0};
  }

private:
  /// Whether the boundary value adds the cell gradient along the tangential part of the
  /// cell-to-face vector
  const bool _two_term_expansion;
  /// The cell gradients of the variable, with the two-term expansion
  Moose::Kokkos::LinearFVCellGradient _gradient;
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVScalarSymmetryBC::BoundaryRelation
KokkosLinearFVScalarSymmetryBC::computeBoundaryValue(const FaceDatum & datum) const
{
  if (!_two_term_expansion)
    return {1, 0};

  const auto normal = datum.faceInfoNormal();
  const auto d_cf = datum.cellToFaceVector();
  const auto correction_vector = d_cf - (d_cf * normal) * normal;
  return {1, _gradient(kokkosSystems(), datum.elem()) * correction_vector};
}
