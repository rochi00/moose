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
 * Provides Dirichlet boundary value and normal-gradient relations for Kokkos linear finite volume
 * kernels. The boundary value is supplied through a Kokkos-compatible functor.
 */
class KokkosLinearFVFunctorDirichletBC : public Moose::Kokkos::LinearFVBoundaryCondition
{
public:
  static InputParameters validParams();

  KokkosLinearFVFunctorDirichletBC(const InputParameters & parameters);

  virtual bool needsBoundaryNonorthogonalCorrection() const override { return true; }

  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryValue(const FaceDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryNormalGradient(const FaceDatum & datum) const;

private:
  /// The functor providing the Dirichlet value on the boundary
  const Moose::Kokkos::ReferenceWrapper<const KokkosParsedFunction> _functor;
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVFunctorDirichletBC::BoundaryRelation
KokkosLinearFVFunctorDirichletBC::computeBoundaryValue(const FaceDatum & datum) const
{
  return {0, _functor->value(_t, datum.faceCentroid())};
}

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVFunctorDirichletBC::BoundaryRelation
KokkosLinearFVFunctorDirichletBC::computeBoundaryNormalGradient(const FaceDatum & datum) const
{
  const auto distance = datum.cellToFaceDistance();
  return {-1 / distance, _functor->value(_t, datum.faceCentroid()) / distance};
}
