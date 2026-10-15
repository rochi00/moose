//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosLinearFVExtrapolatedBC.h"

/**
 * An outflow boundary for Kokkos linear finite volume kernels: the boundary value is extrapolated
 * from the cell as KokkosLinearFVExtrapolatedBC, and there is no diffusive flux. It is the Kokkos
 * analog of LinearFVAdvectionDiffusionOutflowBC.
 */
class KokkosLinearFVOutflowBC : public KokkosLinearFVExtrapolatedBC
{
public:
  static InputParameters validParams();

  KokkosLinearFVOutflowBC(const InputParameters & parameters);

  virtual bool includesMaterialPropertyMultiplier() const override { return true; }

  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryNormalGradient(const FaceDatum &) const
  {
    return {0, 0};
  }
};
