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
 * Provides prescribed outward normal gradient relations for Kokkos linear finite volume kernels.
 */
class KokkosLinearFVFunctorNeumannBC : public Moose::Kokkos::LinearFVBoundaryCondition
{
public:
  static InputParameters validParams();

  KokkosLinearFVFunctorNeumannBC(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryValue(const FaceDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryNormalGradient(const FaceDatum & datum) const;

private:
  /// The functor providing the outward normal gradient on the boundary
  const Moose::Kokkos::ReferenceWrapper<const KokkosParsedFunction> _normal_gradient_functor;
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVFunctorNeumannBC::BoundaryRelation
KokkosLinearFVFunctorNeumannBC::computeBoundaryValue(const FaceDatum & datum) const
{
  return {1, _normal_gradient_functor->value(_t, datum.faceCentroid()) * datum.faceDCFMag()};
}

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVFunctorNeumannBC::BoundaryRelation
KokkosLinearFVFunctorNeumannBC::computeBoundaryNormalGradient(const FaceDatum & datum) const
{
  return {0, _normal_gradient_functor->value(_t, datum.faceCentroid())};
}
