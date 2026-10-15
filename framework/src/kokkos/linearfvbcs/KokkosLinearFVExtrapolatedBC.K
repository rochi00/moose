//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVExtrapolatedBC.h"

registerKokkosLinearFVBoundaryCondition("MooseApp", KokkosLinearFVExtrapolatedBC);

InputParameters
KokkosLinearFVExtrapolatedBC::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVBoundaryCondition::validParams();
  params.addClassDescription(
      "Kokkos LinearFV boundary condition which calculates the face values and face gradients "
      "assuming one or two term expansions from the cell centroid.");
  params.addParam<bool>(
      "use_two_term_expansion",
      false,
      "If an approximate linear expansion should be used to compute the face value.");
  return params;
}

KokkosLinearFVExtrapolatedBC::KokkosLinearFVExtrapolatedBC(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVBoundaryCondition(parameters),
    _two_term_expansion(getParam<bool>("use_two_term_expansion")),
    _u(getLinearFVCellValues(getParam<LinearVariableName>("variable")))
{
  if (_two_term_expansion)
    requestLinearFVCellGradient(getParam<LinearVariableName>("variable"), _gradient);
}
