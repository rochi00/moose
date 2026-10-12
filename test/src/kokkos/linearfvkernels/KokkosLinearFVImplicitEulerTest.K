//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVImplicitEulerTest.h"
#include "FEProblemBase.h"
#include "SystemBase.h"

registerKokkosLinearFVKernel("MooseTestApp", KokkosLinearFVImplicitEulerTest);

InputParameters
KokkosLinearFVImplicitEulerTest::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVElementalKernel::validParams();
  params.addParam<Real>("factor", 1.0, "Constant multiplier on the time derivative.");
  params.addClassDescription("Test-only Kokkos linear FV implicit Euler time derivative.");
  return params;
}

KokkosLinearFVImplicitEulerTest::KokkosLinearFVImplicitEulerTest(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVElementalKernel(parameters),
    _factor(getParam<Real>("factor")),
    _old_tag(
        [this]()
        {
          // Make sure the system keeps the old solution so that its tag exists and is synced
          _fe_problem.getVariable(0, getParam<LinearVariableName>("variable"))
              .sys()
              .needSolutionState(1);
          return _fe_problem.getVectorTagID(Moose::OLD_SOLUTION_TAG);
        }()),
    _dt(TransientInterface::_dt)
{
}

KokkosLinearFVImplicitEulerTest::KokkosLinearFVImplicitEulerTest(
    const KokkosLinearFVImplicitEulerTest & object)
  : Moose::Kokkos::LinearFVElementalKernel(object),
    _factor(object._factor),
    _old_tag(object._old_tag),
    _dt(object._dt)
{
}
