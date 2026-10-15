//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVScalarSymmetryBC.h"

registerKokkosLinearFVBoundaryCondition("MooseApp", KokkosLinearFVScalarSymmetryBC);

InputParameters
KokkosLinearFVScalarSymmetryBC::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVBoundaryCondition::validParams();
  params.addClassDescription("Kokkos LinearFV symmetry boundary condition for a scalar quantity.");
  params.addParam<bool>(
      "use_two_term_expansion",
      false,
      "If an approximate linear expansion should be used to compute the face value.");
  return params;
}

KokkosLinearFVScalarSymmetryBC::KokkosLinearFVScalarSymmetryBC(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVBoundaryCondition(parameters),
    _two_term_expansion(getParam<bool>("use_two_term_expansion"))
{
  if (_two_term_expansion)
    requestLinearFVCellGradient(getParam<LinearVariableName>("variable"), _gradient);
}
