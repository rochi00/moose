//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVLaserBeamFoamLatentHeatAdvection.h"

#include "LinearFVBoundaryCondition.h"
#include "RhieChowMassFlux.h"

registerMooseObject("NavierStokesApp", LinearFVLaserBeamFoamLatentHeatAdvection);
registerMooseObjectAliased("NavierStokesApp",
                           LinearFVLaserBeamFoamLatentHeatAdvection,
                           "LinearFVEnthalpyLiquidFractionLatentHeatAdvection");

InputParameters
LinearFVLaserBeamFoamLatentHeatAdvection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription("Adds the enthalpy liquid-fraction latent heat RHS contribution "
                             "-L * div(rho_face_flux * liquid_fraction).");
  params.addRequiredParam<MooseFunctorName>("L", "Latent heat.");
  params.addRequiredParam<MooseFunctorName>("liquid_fraction", "Liquid fraction.");
  params.addParam<UserObjectName>(
      "rhie_chow_user_object", "", "The Rhie-Chow user object providing the mass face flux.");
  params.addParam<MooseFunctorName>(
      "mass_flux_functor",
      "",
      "Optional face mass-flux functor. When supplied, this kernel uses it directly instead of "
      "querying the Rhie-Chow user object.");
  params.addParam<bool>(
      "mass_flux_is_integrated",
      false,
      "Whether mass_flux_functor already includes the face area. When true, the kernel divides by "
      "the current face area before applying the usual linear FV face-area scaling.");
  return params;
}

LinearFVLaserBeamFoamLatentHeatAdvection::LinearFVLaserBeamFoamLatentHeatAdvection(
    const InputParameters & params)
  : LinearFVFluxKernel(params),
    _L(getFunctor<Real>("L")),
    _liquid_fraction(getFunctor<Real>("liquid_fraction")),
    _mass_flux_provider(isParamValid("rhie_chow_user_object") &&
                                !getParam<UserObjectName>("rhie_chow_user_object").empty()
                            ? &getUserObject<RhieChowMassFlux>("rhie_chow_user_object")
                            : nullptr),
    _mass_flux_functor(isParamValid("mass_flux_functor") &&
                               !getParam<MooseFunctorName>("mass_flux_functor").empty()
                           ? &getFunctor<Real>("mass_flux_functor")
                           : nullptr),
    _mass_flux_is_integrated(getParam<bool>("mass_flux_is_integrated"))
{
  if (!_mass_flux_provider && !_mass_flux_functor)
    paramError("rhie_chow_user_object",
               "Either 'rhie_chow_user_object' or 'mass_flux_functor' must be provided.");
}

Real
LinearFVLaserBeamFoamLatentHeatAdvection::computeElemMatrixContribution()
{
  return 0.0;
}

Real
LinearFVLaserBeamFoamLatentHeatAdvection::computeNeighborMatrixContribution()
{
  return 0.0;
}

Real
LinearFVLaserBeamFoamLatentHeatAdvection::computeElemRightHandSideContribution()
{
  return -computeFaceContribution();
}

Real
LinearFVLaserBeamFoamLatentHeatAdvection::computeNeighborRightHandSideContribution()
{
  return computeFaceContribution();
}

Real
LinearFVLaserBeamFoamLatentHeatAdvection::computeBoundaryMatrixContribution(
    const LinearFVBoundaryCondition & /*bc*/)
{
  return 0.0;
}

Real
LinearFVLaserBeamFoamLatentHeatAdvection::computeBoundaryRHSContribution(
    const LinearFVBoundaryCondition & /*bc*/)
{
  return -computeFaceContribution();
}

Real
LinearFVLaserBeamFoamLatentHeatAdvection::computeFaceContribution()
{
  const auto state = determineState();

  if (!_cached_rhs_contribution)
  {
    _cached_rhs_contribution = true;

    const auto face_arg = functorFaceArg(_liquid_fraction, *_current_face_info);
    const Real latent_heat = _L(face_arg, state);
    const Real liquid_fraction = _liquid_fraction(face_arg, state);
    Real rho_flux = 0.0;
    if (_mass_flux_functor)
    {
      rho_flux = (*_mass_flux_functor)(face_arg, state);
      if (_mass_flux_is_integrated)
        rho_flux = _current_face_area > 0.0 ? rho_flux / _current_face_area : 0.0;
    }
    else
      rho_flux = _mass_flux_provider->getMassFlux(*_current_face_info);

    _face_contribution = latent_heat * rho_flux * liquid_fraction * _current_face_area;
  }

  return _face_contribution;
}
