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
 * Test-only Kokkos linear FV source equal to another linear FV variable, which may belong to
 * another system: rhs += coupled * V. Reproduces LinearFVSource with source_density set to the
 * variable.
 */
class KokkosLinearFVCoupledSourceTest : public Moose::Kokkos::LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  KokkosLinearFVCoupledSourceTest(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION CellContribution computeCellContribution(const CellDatum & datum) const
  {
    return {0, _coupled(kokkosSystems(), datum.elem()) * datum.elemVolume()};
  }

private:
  /// The cell values of the coupled variable
  const Moose::Kokkos::LinearFVCellValues _coupled;
};
