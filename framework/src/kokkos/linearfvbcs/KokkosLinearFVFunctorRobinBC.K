//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVFunctorRobinBC.h"

#include "MooseUtils.h"

registerKokkosLinearFVBoundaryCondition("MooseApp", KokkosLinearFVFunctorRobinBC);

InputParameters
KokkosLinearFVFunctorRobinBC::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVBoundaryCondition::validParams();
  params.addClassDescription(
      "Kokkos LinearFV Robin boundary condition of the form alpha * grad(u).n + beta * u = gamma "
      "with functor coefficients.");
  params.addParam<MooseFunctorName>(
      "alpha", 1.0, "Functor for the coefficient of the normal gradient term.");
  params.addParam<MooseFunctorName>("beta", 1.0, "Functor for the coefficient of the scalar term.");
  params.addParam<MooseFunctorName>(
      "gamma", 1.0, "Functor for the constant term on the RHS of the Robin BC.");
  return params;
}

KokkosLinearFVFunctorRobinBC::KokkosLinearFVFunctorRobinBC(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVBoundaryCondition(parameters)
{
  Real alpha;
  if (MooseUtils::parsesToReal(getParam<MooseFunctorName>("alpha"), &alpha) &&
      MooseUtils::isZero(alpha))
    paramError("alpha", "This value shall not be 0. Use a Dirichlet boundary condition instead!");

  // As the host, the coefficients are evaluated from the variable's side of the face
  PropertyUse use;
  use.boundary_faces = true;
  requestLinearFVProperty("alpha", _alpha, use);
  requestLinearFVProperty("beta", _beta, use);
  requestLinearFVProperty("gamma", _gamma, use);
  requestLinearFVCellGradient(getParam<LinearVariableName>("variable"), _gradient);
}
