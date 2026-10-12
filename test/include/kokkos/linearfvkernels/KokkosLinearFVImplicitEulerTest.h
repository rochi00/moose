//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosLinearFVElementalKernel.h"

/**
 * Test-only Kokkos linear FV implicit Euler time derivative with a constant factor:
 *   matrix += factor * V / dt,  rhs += factor * V / dt * u_old
 * Reproduces LinearFVTimeDerivative with ImplicitEuler and a constant factor. It reads the old
 * solution through the OLD_SOLUTION_TAG vector of the Kokkos system.
 */
class KokkosLinearFVImplicitEulerTest : public Moose::Kokkos::LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  KokkosLinearFVImplicitEulerTest(const InputParameters & parameters);
  KokkosLinearFVImplicitEulerTest(const KokkosLinearFVImplicitEulerTest & object);

  template <typename Derived>
  KOKKOS_FUNCTION Real computeMatrixContribution(const Moose::Kokkos::FVDatum & datum) const;

  template <typename Derived>
  KOKKOS_FUNCTION Real computeRightHandSideContribution(const Moose::Kokkos::FVDatum & datum) const;

private:
  /// Constant multiplier on the time derivative
  const Real _factor;
  /// Vector tag of the old solution
  const TagID _old_tag;
  /// Time step size, refreshed at every dispatch
  Moose::Kokkos::Scalar<Real> _dt;
};

template <typename Derived>
KOKKOS_FUNCTION inline Real
KokkosLinearFVImplicitEulerTest::computeMatrixContribution(
    const Moose::Kokkos::FVDatum & datum) const
{
  return _factor * datum.elementVolume() / _dt;
}

template <typename Derived>
KOKKOS_FUNCTION inline Real
KokkosLinearFVImplicitEulerTest::computeRightHandSideContribution(
    const Moose::Kokkos::FVDatum & datum) const
{
  const auto & sys = kokkosSystem(_var.sys());
  const auto dof = sys.getElemLocalDofIndex(datum.elemID(), 0, _var.var());
  return _factor * datum.elementVolume() / _dt * sys.getVectorDofValue(dof, _old_tag);
}
