//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CircleAreaFractionFunction.h"

registerMooseObject("NavierStokesApp", CircleAreaFractionFunction);

InputParameters
CircleAreaFractionFunction::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription(
      "Computes the area fraction of each cell occupied by a circle, "
      "using subgrid sampling. Equivalent to Basilisk's fraction() for circles.");
  params.addRequiredParam<Real>("radius", "Circle radius.");
  params.addRequiredParam<Real>("x_center", "Circle center x-coordinate.");
  params.addRequiredParam<Real>("y_center", "Circle center y-coordinate.");
  params.addRequiredParam<Real>("cell_size", "Cell size dx (uniform mesh).");
  params.addParam<unsigned int>(
      "n_subdivisions", 50, "Number of sub-cells per direction for sampling.");
  return params;
}

CircleAreaFractionFunction::CircleAreaFractionFunction(const InputParameters & params)
  : Function(params),
    _radius(getParam<Real>("radius")),
    _x_center(getParam<Real>("x_center")),
    _y_center(getParam<Real>("y_center")),
    _dx(getParam<Real>("cell_size")),
    _n_subdivisions(getParam<unsigned int>("n_subdivisions"))
{
}

Real
CircleAreaFractionFunction::value(Real /*t*/, const Point & p) const
{
  const Real R_sq = _radius * _radius;

  // Quick check: if cell is entirely inside or outside the circle, skip sampling
  const Real dist = std::sqrt((p(0) - _x_center) * (p(0) - _x_center) +
                              (p(1) - _y_center) * (p(1) - _y_center));
  // Cell diagonal half-length
  const Real half_diag = _dx * std::sqrt(2.0) / 2.0;

  if (dist + half_diag <= _radius)
    return 1.0; // entirely inside
  if (dist - half_diag >= _radius)
    return 0.0; // entirely outside

  // Cell lower-left corner
  const Real x0 = p(0) - _dx / 2.0;
  const Real y0 = p(1) - _dx / 2.0;
  const Real sub_dx = _dx / _n_subdivisions;

  unsigned int count = 0;
  for (unsigned int i = 0; i < _n_subdivisions; ++i)
  {
    const Real xs = x0 + (i + 0.5) * sub_dx;
    const Real dxs = xs - _x_center;
    for (unsigned int j = 0; j < _n_subdivisions; ++j)
    {
      const Real ys = y0 + (j + 0.5) * sub_dx;
      const Real dys = ys - _y_center;
      if (dxs * dxs + dys * dys <= R_sq)
        count++;
    }
  }

  return static_cast<Real>(count) / static_cast<Real>(_n_subdivisions * _n_subdivisions);
}
