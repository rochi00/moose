//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVLevelSetAdvection.h"
#include "MooseLinearVariableFV.h"
#include "NSFVUtils.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", LinearFVLevelSetAdvection);

InputParameters
LinearFVLevelSetAdvection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Advection kernel for the level-set variable using the volumetric face flux "
      "from RhieChowMassFluxMultiPhase. Standard upwind/TVD interpolation without "
      "MULES or compression velocity.");
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object",
      "The rhie-chow user-object which provides the volumetric face flux.");
  params += Moose::FV::advectedInterpolationParameter();
  return params;
}

LinearFVLevelSetAdvection::LinearFVLevelSetAdvection(const InputParameters & params)
  : LinearFVFluxKernel(params),
    _mass_flux_provider(getUserObject<RhieChowMassFluxMultiPhase>("rhie_chow_user_object")),
    _advected_interp_coeffs(std::make_pair<Real, Real>(0, 0)),
    _volumetric_face_flux(0.0)
{
  Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");
}

Real
LinearFVLevelSetAdvection::computeElemMatrixContribution()
{
  return _advected_interp_coeffs.first * _volumetric_face_flux * _current_face_area;
}

Real
LinearFVLevelSetAdvection::computeNeighborMatrixContribution()
{
  return _advected_interp_coeffs.second * _volumetric_face_flux * _current_face_area;
}

Real
LinearFVLevelSetAdvection::computeElemRightHandSideContribution()
{
  return 0.0;
}

Real
LinearFVLevelSetAdvection::computeNeighborRightHandSideContribution()
{
  return 0.0;
}

Real
LinearFVLevelSetAdvection::computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  const auto boundary_value_matrix_contrib = adv_bc->computeBoundaryValueMatrixContribution();
  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  return boundary_value_matrix_contrib * factor * _volumetric_face_flux * _current_face_area;
}

Real
LinearFVLevelSetAdvection::computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0);
  const auto boundary_value_rhs_contrib = adv_bc->computeBoundaryValueRHSContribution();

  return -boundary_value_rhs_contrib * factor * _volumetric_face_flux * _current_face_area;
}

void
LinearFVLevelSetAdvection::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  _volumetric_face_flux = _mass_flux_provider.getVolumetricFaceFlux(*face_info);

  _advected_interp_coeffs =
      interpCoeffs(_advected_interp_method, *_current_face_info, true, _volumetric_face_flux);
}
