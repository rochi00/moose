//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "MooseTypes.h"

#include <algorithm>
#include <cmath>

namespace LowMachCUI
{
/**
 * Koren's convection-bounded cubic-upwind face reconstruction.
 *
 * The arguments follow the normalized-variable convention: U is far-upwind, C is upwind, and D
 * is downwind.
 */
inline Real
faceValue(const Real upwind_value, const Real far_upwind_value, const Real downwind_value)
{
  const Real denominator = downwind_value - far_upwind_value;
  const Real scale =
      std::max({1.0, std::abs(upwind_value), std::abs(far_upwind_value), std::abs(downwind_value)});
  if (std::abs(denominator) <= TOLERANCE * scale)
    return upwind_value;

  const Real normalized_upwind = (upwind_value - far_upwind_value) / denominator;
  Real normalized_face = normalized_upwind;
  if (normalized_upwind > 0.0 && normalized_upwind <= 2.0 / 13.0)
    normalized_face = 3.0 * normalized_upwind;
  else if (normalized_upwind > 2.0 / 13.0 && normalized_upwind <= 4.0 / 5.0)
    normalized_face = 5.0 / 6.0 * normalized_upwind + 1.0 / 3.0;
  else if (normalized_upwind > 4.0 / 5.0 && normalized_upwind <= 1.0)
    normalized_face = 1.0;

  return far_upwind_value + normalized_face * denominator;
}
}
