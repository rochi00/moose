//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearWCNSFV2PMassDriftFlux.h"
#include "NS.h"
#include "NavierStokesMethods.h"

registerMooseObject("NavierStokesApp", LinearWCNSFV2PMassDriftFlux);

InputParameters
LinearWCNSFV2PMassDriftFlux::validParams()
{
  InputParameters params = LinearWCNSFV2PDriftFluxBase::validParams();
  params.addClassDescription(
      "Adds the dilatation produced by the relative motion of the phases to the pressure equation "
      "of a two phase mixture, so that the mass averaged mixture velocity satisfies its own "
      "continuity equation rather than being constrained to be solenoidal.");
  params.makeParamRequired<MooseFunctorName>("fraction_dispersed");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Density of the mixture.");
  return params;
}

LinearWCNSFV2PMassDriftFlux::LinearWCNSFV2PMassDriftFlux(const InputParameters & params)
  : LinearWCNSFV2PDriftFluxBase(params),
    _rho_mixture(getFunctor<Real>(NS::density)),
    _face_drift_flux(0.0)
{
}

Real
LinearWCNSFV2PMassDriftFlux::computeDriftVolumetricFlux()
{
  const auto * const fi = _current_face_info;

  // The phases cannot separate across an impermeable boundary, so no dilatation is produced there
  if (!slipAllowedOnCurrentFace())
    return 0.0;

  const auto state = determineState();
  const auto face_arg = makeCDFace(*fi);

  const auto alpha = _f_d(face_arg, state);
  const auto rho_m = _rho_mixture(face_arg, state);

  // c_d is the mass fraction of the dispersed phase. The factor alpha - c_d vanishes when the two
  // phase densities are equal, which is what makes this term disappear in a single density mixture.
  const auto c_d = (rho_m > 0.0) ? alpha * _rho_d(face_arg, state) / rho_m : 0.0;
  const auto factor = alpha - c_d;

  const auto slip = slipVelocity(face_arg, state);
  return factor * (slip * fi->normal());
}

void
LinearWCNSFV2PMassDriftFlux::setupFaceData(const FaceInfo * face_info)
{
  LinearWCNSFV2PDriftFluxBase::setupFaceData(face_info);
  _face_drift_flux = computeDriftVolumetricFlux();
}

Real
LinearWCNSFV2PMassDriftFlux::computeElemMatrixContribution()
{
  return 0.0;
}

Real
LinearWCNSFV2PMassDriftFlux::computeNeighborMatrixContribution()
{
  return 0.0;
}

Real
LinearWCNSFV2PMassDriftFlux::computeElemRightHandSideContribution()
{
  return _face_drift_flux * _current_face_area;
}

Real
LinearWCNSFV2PMassDriftFlux::computeNeighborRightHandSideContribution()
{
  return -_face_drift_flux * _current_face_area;
}

Real
LinearWCNSFV2PMassDriftFlux::computeBoundaryMatrixContribution(
    const LinearFVBoundaryCondition & /*bc*/)
{
  return 0.0;
}

Real
LinearWCNSFV2PMassDriftFlux::computeBoundaryRHSContribution(
    const LinearFVBoundaryCondition & /*bc*/)
{
  return _face_drift_flux * _current_face_area;
}
