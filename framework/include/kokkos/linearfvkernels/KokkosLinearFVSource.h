//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosLinearFVElementalKernel.h"

/**
 * Kokkos linear finite volume elemental kernel that adds a scaled volumetric source to the
 * right-hand side. It is the Kokkos analog of LinearFVSource.
 */
class KokkosLinearFVSource : public Moose::Kokkos::LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  KokkosLinearFVSource(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION CellContribution computeCellContribution(const CellDatum & datum) const;

private:
  /// The source density
  Moose::Kokkos::LinearFVProperty _source_density;
  /// Coefficient for multiplying the source density
  Moose::Kokkos::LinearFVProperty _scale;
};

template <typename Derived>
KOKKOS_FUNCTION inline KokkosLinearFVSource::CellContribution
KokkosLinearFVSource::computeCellContribution(const CellDatum & datum) const
{
  const auto elem = datum.elem();
  return {0, _scale.cell(elem) * _source_density.cell(elem) * datum.elemVolume()};
}
