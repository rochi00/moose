//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVDiffusion.h"

#include "FEProblemBase.h"

registerKokkosLinearFVKernel("MooseApp", KokkosLinearFVDiffusion);

InputParameters
KokkosLinearFVDiffusion::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVFluxKernel::validParams();
  params.addParam<MooseFunctorName>("diffusion_coeff", 1.0, "The diffusion coefficient.");
  params.addParam<InterpolationMethodName>(
      "coeff_interp_method",
      "Optional finite volume interpolation method used to compute a face-centered diffusion "
      "coefficient. If omitted, the functor is evaluated directly on the face.");
  params.addParam<bool>(
      "use_nonorthogonal_correction",
      true,
      "If the nonorthogonal correction should be used when computing the normal gradient.");
  params.addClassDescription("Kokkos LinearFV diffusion kernel.");
  return params;
}

KokkosLinearFVDiffusion::KokkosLinearFVDiffusion(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVFluxKernel(parameters),
    _use_nonorthogonal_correction(getParam<bool>("use_nonorthogonal_correction"))
{
  PropertyUse use;
  use.internal_faces = true;
  use.boundary_faces = true;
  if (isParamValid("coeff_interp_method"))
    use.interpolation = &_fe_problem.getFVFaceInterpolationMethod(
        getParam<InterpolationMethodName>("coeff_interp_method"));
  requestLinearFVProperty("diffusion_coeff", _diffusion_coeff, use);

  if (_use_nonorthogonal_correction)
    requestLinearFVCellGradient(getParam<LinearVariableName>("variable"), _gradient);
}
