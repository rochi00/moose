//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosFVGreenGaussGradient.h"
#include "KokkosGreenGaussGradient.h"

registerMooseObject("MooseApp", KokkosFVGreenGaussGradient);

InputParameters
KokkosFVGreenGaussGradient::validParams()
{
  InputParameters params = FVGradientMethod::validParams();
  params.addClassDescription("Green-Gauss cell-centered gradient method computed with Kokkos.");
  return params;
}

KokkosFVGreenGaussGradient::KokkosFVGreenGaussGradient(const InputParameters & params)
  : FVGradientMethod(params), _gradient(std::make_unique<Moose::Kokkos::GreenGaussGradient>())
{
}

KokkosFVGreenGaussGradient::~KokkosFVGreenGaussGradient() = default;

void
KokkosFVGreenGaussGradient::computeGradientWithoutLimiter(
    SystemBase & system,
    GradientContainer & gradient,
    const std::unordered_set<unsigned int> & variable_numbers) const
{
  _gradient->compute(system, gradient, variable_numbers);
}

void
KokkosFVGreenGaussGradient::applyLimiter(
    SystemBase & system,
    GradientContainer & gradient,
    const std::unordered_set<unsigned int> & variable_numbers) const
{
  _gradient->limit(system, gradient, variable_numbers, limiterType());
}
