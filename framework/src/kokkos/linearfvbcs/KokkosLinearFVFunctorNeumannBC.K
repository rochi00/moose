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
  params.addRequiredParam<MooseFunctorName>(
      "functor", "The outward diffusive flux functor on the boundary.");
  params.addParam<MooseFunctorName>("diffusion_coeff", 1.0, "The diffusion coefficient.");
  params.addClassDescription("Kokkos LinearFV Neumann boundary condition with a prescribed "
                             "diffusive flux.");
  return params;
}

KokkosLinearFVFunctorNeumannBC::KokkosLinearFVFunctorNeumannBC(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVBoundaryCondition(parameters)
{
  PropertyUse use;
  use.boundary_faces = true;
  requestLinearFVProperty("functor", _flux_functor, use);
  requestLinearFVProperty("diffusion_coeff", _diffusion_coeff, use);
  requestLinearFVCellGradient(getParam<LinearVariableName>("variable"), _gradient);
}
