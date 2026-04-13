//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Function.h"

/**
 * Computes the area fraction of each cell occupied by a circle.
 * Equivalent to Basilisk's fraction() function for a circle.
 *
 * Uses subgrid sampling: divides each cell into NxN sub-cells and counts
 * the fraction inside the circle. With N=50, accuracy is ~4e-4 per cell,
 * which gives O(dx^2) curvature convergence with height functions.
 *
 * Assumes a uniform Cartesian mesh with square cells of size dx.
 */
class CircleAreaFractionFunction : public Function
{
public:
  static InputParameters validParams();
  CircleAreaFractionFunction(const InputParameters & params);

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

protected:
  /// Circle radius
  const Real _radius;
  /// Circle center x
  const Real _x_center;
  /// Circle center y
  const Real _y_center;
  /// Cell size (uniform mesh)
  const Real _dx;
  /// Number of sub-cells per direction for sampling
  const unsigned int _n_subdivisions;
};
