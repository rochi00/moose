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
 * Extrapolates the boundary value from the cell value, optionally with the cell gradient, for
 * Kokkos linear finite volume kernels. It is the Kokkos analog of
 * LinearFVAdvectionDiffusionExtrapolatedBC.
 */
class KokkosLinearFVExtrapolatedBC : public Moose::Kokkos::LinearFVBoundaryCondition
{
public:
  static InputParameters validParams();

  KokkosLinearFVExtrapolatedBC(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryValue(const FaceDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryNormalGradient(const FaceDatum & datum) const;

protected:
  /// The gradient term of the boundary value: the cell gradient along the cell-to-face vector with
  /// the two-term expansion, else zero
  KOKKOS_FUNCTION Real gradientTerm(const FaceDatum & datum) const;

  /// Whether the boundary value adds the cell gradient along the cell-to-face vector
  const bool _two_term_expansion;
  /// The cell values of the variable
  const Moose::Kokkos::LinearFVCellValues _u;
  /// The cell gradients of the variable, with the two-term expansion
  Moose::Kokkos::LinearFVCellGradient _gradient;
};

KOKKOS_FUNCTION inline Real
KokkosLinearFVExtrapolatedBC::gradientTerm(const FaceDatum & datum) const
{
  return _two_term_expansion ? _gradient(kokkosSystems(), datum.elem()) * datum.cellToFaceVector()
                             : 0;
}

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVExtrapolatedBC::BoundaryRelation
KokkosLinearFVExtrapolatedBC::computeBoundaryValue(const FaceDatum & datum) const
{
  return {1, gradientTerm(datum)};
}

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVExtrapolatedBC::BoundaryRelation
KokkosLinearFVExtrapolatedBC::computeBoundaryNormalGradient(const FaceDatum & datum) const
{
  // As the host, the difference to the boundary value computed with the current cell value, which
  // the right-hand side carries
  const auto distance = datum.cellToFaceDistance();
  const auto boundary_value = _u(kokkosSystems(), datum.elem()) + gradientTerm(datum);
  return {-1 / distance, boundary_value / distance};
}
