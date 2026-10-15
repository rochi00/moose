//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVFunctorNeumannBC.h"

registerKokkosLinearFVBoundaryCondition("MooseApp", KokkosLinearFVFunctorNeumannBC);

InputParameters
KokkosLinearFVFunctorNeumannBC::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVBoundaryCondition::validParams();
  params.addRequiredParam<FunctionName>("functor",
                                        "The outward diffusive flux functor on the boundary.");
  params.addRequiredParam<FunctionName>("diffusion_coeff", "The diffusion coefficient.");
  params.addClassDescription("Kokkos LinearFV Neumann boundary condition with a prescribed "
                             "diffusive flux.");
  return params;
}

KokkosLinearFVFunctorNeumannBC::KokkosLinearFVFunctorNeumannBC(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVBoundaryCondition(parameters),
    _flux_functor(getKokkosFunction<KokkosParsedFunction>("functor")),
    _diffusion_coeff(getKokkosFunction<KokkosParsedFunction>("diffusion_coeff"))
{
  requestLinearFVCellGradient(getParam<LinearVariableName>("variable"), _gradient);
}
