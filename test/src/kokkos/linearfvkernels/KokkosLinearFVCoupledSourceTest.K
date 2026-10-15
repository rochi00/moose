//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVCoupledSourceTest.h"

registerKokkosLinearFVKernel("MooseTestApp", KokkosLinearFVCoupledSourceTest);

InputParameters
KokkosLinearFVCoupledSourceTest::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVElementalKernel::validParams();
  params.addRequiredParam<VariableName>("coupled", "The linear FV variable used as the source.");
  params.addClassDescription("Test-only Kokkos linear FV source equal to a linear FV variable.");
  return params;
}

KokkosLinearFVCoupledSourceTest::KokkosLinearFVCoupledSourceTest(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVElementalKernel(parameters),
    _coupled(getLinearFVCellValues(getParam<VariableName>("coupled")))
{
}
