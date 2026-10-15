//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVAdvection.h"

#include "FEProblemBase.h"
#include "FVAdvectedInterpolationMethod.h"
#include "FVAdvectedMinmodWeightBased.h"
#include "FVAdvectedMUSCLDeferredCorrection.h"
#include "FVAdvectedUpwind.h"
#include "FVAdvectedVanLeerWeightBased.h"
#include "FVGeometricAverage.h"

registerKokkosLinearFVKernel("MooseApp", KokkosLinearFVAdvection);

InputParameters
KokkosLinearFVAdvection::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVFluxKernel::validParams();
  params.addRequiredParam<RealVectorValue>("velocity", "Constant advection velocity.");
  params.addRequiredParam<InterpolationMethodName>(
      "advected_interp_method_name",
      "Name of the FVInterpolationMethod to use for the advected quantity.");
  params.addClassDescription("Kokkos LinearFV advection kernel.");
  return params;
}

KokkosLinearFVAdvection::KokkosLinearFVAdvection(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVFluxKernel(parameters),
    _velocity(getParam<RealVectorValue>("velocity")),
    _u(getLinearFVCellValues(getParam<LinearVariableName>("variable")))
{
  const auto & method = _fe_problem.getFVAdvectedInterpolationMethod(
      getParam<InterpolationMethodName>("advected_interp_method_name"), 0);
  const auto & method_object = dynamic_cast<const MooseObject &>(method);

  if (dynamic_cast<const FVAdvectedUpwind *>(&method))
    _interpolation = Interpolation::UPWIND;
  else if (dynamic_cast<const FVGeometricAverage *>(&method))
    _interpolation = Interpolation::AVERAGE;
  else if (dynamic_cast<const FVAdvectedMinmodWeightBased *>(&method))
    _interpolation = Interpolation::MINMOD;
  else if (dynamic_cast<const FVAdvectedVanLeerWeightBased *>(&method))
    _interpolation = Interpolation::VAN_LEER;
  else if (dynamic_cast<const FVAdvectedMUSCLDeferredCorrection *>(&method))
    _interpolation = Interpolation::MUSCL;
  else
    paramError("advected_interp_method_name",
               "The interpolation method type '",
               method_object.type(),
               "' is not supported by ",
               type(),
               ".");

  if (_interpolation == Interpolation::MINMOD || _interpolation == Interpolation::VAN_LEER)
  {
    _limit_to_linear = method_object.getParam<bool>("limit_to_linear");
    _blending_factor = method_object.getParam<Real>("blending_factor");
  }
  else if (_interpolation == Interpolation::MUSCL)
    _deferred_correction_factor = method_object.getParam<Real>("deferred_correction_factor");

  if (method.needsGradients())
    requestLinearFVCellGradient(
        getParam<LinearVariableName>("variable"), _gradient, method.gradientMethodName());
}
