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

registerMooseObject("NavierStokesApp", LinearWCNSFV2PMassDriftFlux);

InputParameters
LinearWCNSFV2PMassDriftFlux::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Adds the dilatation produced by the relative motion of the phases to the pressure equation "
      "of a two phase mixture, so that the mass averaged mixture velocity satisfies its own "
      "continuity equation rather than being constrained to be solenoidal.");
  params.addRequiredParam<MooseFunctorName>("fraction_dispersed",
                                            "Volume fraction of the dispersed phase.");
  params.addRequiredParam<MooseFunctorName>("rho_d", "Density of the dispersed phase.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Density of the mixture.");
  params.addRequiredParam<MooseFunctorName>("u_slip", "The slip velocity in the x direction.");
  params.addParam<MooseFunctorName>("v_slip", "The slip velocity in the y direction.");
  params.addParam<MooseFunctorName>("w_slip", "The slip velocity in the z direction.");
  params.addParam<std::vector<BoundaryName>>(
      "slip_boundaries",
      {},
      "Boundaries across which the dispersed phase may travel. The phases cannot separate across "
      "an impermeable wall, so the relative motion produces no dilatation there.");
  return params;
}

LinearWCNSFV2PMassDriftFlux::LinearWCNSFV2PMassDriftFlux(const InputParameters & params)
  : LinearFVFluxKernel(params),
    _f_d(getFunctor<Real>("fraction_dispersed")),
    _rho_d(getFunctor<Real>("rho_d")),
    _rho_mixture(getFunctor<Real>(NS::density)),
    _u_slip(getFunctor<Real>("u_slip")),
    _v_slip(isParamValid("v_slip") ? &getFunctor<Real>("v_slip") : nullptr),
    _w_slip(isParamValid("w_slip") ? &getFunctor<Real>("w_slip") : nullptr),
    _face_drift_flux(0.0)
{
  for (const auto & bname : getParam<std::vector<BoundaryName>>("slip_boundaries"))
    _slip_boundaries.insert(_mesh.getBoundaryID(bname));
}

Real
LinearWCNSFV2PMassDriftFlux::computeDriftVolumetricFlux()
{
  const auto * const fi = _current_face_info;

  // The phases cannot separate across an impermeable boundary, so no dilatation is produced there
  const bool slip_allowed_here =
      fi->neighborPtr() ||
      (!fi->boundaryIDs().empty() && _slip_boundaries.count(*fi->boundaryIDs().begin()));

  if (!slip_allowed_here)
    return 0.0;

  const auto state = determineState();
  const auto face_arg = makeCDFace(*fi);

  const auto alpha = _f_d(face_arg, state);
  const auto rho_m = _rho_mixture(face_arg, state);

  // c_d is the mass fraction of the dispersed phase. The factor alpha - c_d vanishes when the two
  // phase densities are equal, which is what makes this term disappear in a single density mixture.
  const auto c_d = (rho_m > 0.0) ? alpha * _rho_d(face_arg, state) / rho_m : 0.0;
  const auto factor = alpha - c_d;

  RealVectorValue slip;
  slip(0) = _u_slip(face_arg, state);
  if (_v_slip)
    slip(1) = (*_v_slip)(face_arg, state);
  if (_w_slip)
    slip(2) = (*_w_slip)(face_arg, state);

  return factor * (slip * fi->normal());
}

void
LinearWCNSFV2PMassDriftFlux::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);
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
