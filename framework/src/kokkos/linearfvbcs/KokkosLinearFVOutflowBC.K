//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVOutflowBC.h"

registerKokkosLinearFVBoundaryCondition("MooseApp", KokkosLinearFVOutflowBC);

InputParameters
KokkosLinearFVOutflowBC::validParams()
{
  InputParameters params = KokkosLinearFVExtrapolatedBC::validParams();
  params.addClassDescription("Kokkos LinearFV boundary condition which represents a surface with "
                             "outflowing material with a constant velocity.");
  return params;
}

KokkosLinearFVOutflowBC::KokkosLinearFVOutflowBC(const InputParameters & parameters)
  : KokkosLinearFVExtrapolatedBC(parameters)
{
}
