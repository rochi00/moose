//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVTimeDerivative.h"

#include "BDF2.h"
#include "ImplicitEuler.h"
#include "SystemBase.h"

registerKokkosLinearFVKernel("MooseApp", KokkosLinearFVTimeDerivative);

InputParameters
KokkosLinearFVTimeDerivative::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVElementalKernel::validParams();
  params.addClassDescription("Kokkos LinearFV time derivative with a constant multiplier.");
  params.addParam<Real>("factor", 1.0, "A constant multiplier on the time derivative.");
  return params;
}

KokkosLinearFVTimeDerivative::KokkosLinearFVTimeDerivative(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVElementalKernel(parameters),
    _factor(getParam<Real>("factor")),
    _time_integrator(getCheckedPointerParam<SystemBase *>("_sys")->getTimeIntegrator(_var.var())),
    _u_old(getLinearFVCellValues(getParam<LinearVariableName>("variable"), 1)),
    _u_older(dynamic_cast<const BDF2 *>(&_time_integrator)
                 ? getLinearFVCellValues(getParam<LinearVariableName>("variable"), 2)
                 : Moose::Kokkos::LinearFVCellValues{})
{
  if (!dynamic_cast<const ImplicitEuler *>(&_time_integrator) &&
      !dynamic_cast<const BDF2 *>(&_time_integrator))
    mooseError("Kokkos linear FV time derivatives support the ImplicitEuler and BDF2 time "
               "integrators; the variable uses ",
               _time_integrator.type(),
               ".");
}

void
KokkosLinearFVTimeDerivative::dispatchContributions()
{
  // The same coefficients as ImplicitEuler and BDF2::timeDerivativeMatrixContribution() and
  // timeDerivativeRHSContribution(); BDF2 takes an implicit Euler first step
  const auto * const bdf2 = dynamic_cast<const BDF2 *>(&_time_integrator);
  if (bdf2 && _t_step > 1)
  {
    const auto & weights = bdf2->weights();
    _matrix_coefficient = weights[0] / _dt;
    _old_coefficient = -weights[1] / _dt;
    _older_coefficient = -weights[2] / _dt;
  }
  else
  {
    _matrix_coefficient = 1 / _dt;
    _old_coefficient = 1 / _dt;
    _older_coefficient = 0;
  }

  Moose::Kokkos::LinearFVElementalKernel::dispatchContributions();
}
