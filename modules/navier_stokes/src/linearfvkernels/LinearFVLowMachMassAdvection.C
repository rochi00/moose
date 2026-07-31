//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LinearFVLowMachMassAdvection.h"

#include "LinearFVAdvectionDiffusionBC.h"
#include "NSFVUtils.h"
#include "RhieChowMassFlux.h"

#include <cmath>

registerMooseObject("NavierStokesApp", LinearFVLowMachMassAdvection);

InputParameters
LinearFVLowMachMassAdvection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Implicitly advects temporary low-Mach density with bounded donor-cell upwinding.");
  params.addParam<MooseFunctorName>(
      "volumetric_face_flux",
      "",
      "Volumetric face-flux density, oriented with the FaceInfo normal.");
  params.addParam<UserObjectName>(
      "rhie_chow_user_object", "", "Rhie-Chow user object providing the volumetric face flux.");
  params.addParam<MooseFunctorName>(
      "material_compression_flux",
      "Optional limited artificial material-compression face-flux density.");
  params.addParam<MooseFunctorName>(
      "material_density",
      "PCM/material-phase density carried by the optional material-compression flux.");
  params.addParam<MooseFunctorName>(
      "background_density",
      "Gas/background-phase density displaced by the optional material-compression flux.");
  return params;
}

LinearFVLowMachMassAdvection::LinearFVLowMachMassAdvection(const InputParameters & params)
  : LinearFVFluxKernel(params),
    _volumetric_face_flux_name(getParam<MooseFunctorName>("volumetric_face_flux")),
    _volumetric_flux_provider(isParamValid("rhie_chow_user_object") &&
                                      !getParam<UserObjectName>("rhie_chow_user_object").empty()
                                  ? &getUserObject<RhieChowMassFlux>("rhie_chow_user_object")
                                  : nullptr),
    _volumetric_face_flux(isParamValid("volumetric_face_flux") &&
                                  !getParam<MooseFunctorName>("volumetric_face_flux").empty()
                              ? &getFunctor<Real>("volumetric_face_flux")
                              : nullptr),
    _material_compression_flux(isParamValid("material_compression_flux")
                                   ? &getFunctor<Real>("material_compression_flux")
                                   : nullptr),
    _material_density(isParamValid("material_density") ? &getFunctor<Real>("material_density")
                                                       : nullptr),
    _background_density(isParamValid("background_density") ? &getFunctor<Real>("background_density")
                                                           : nullptr),
    _face_volumetric_flux(0.0),
    _density_interp_coeffs(0.0, 0.0),
    _face_explicit_correction(0.0)
{
  if (!_volumetric_flux_provider && !_volumetric_face_flux)
    paramError("volumetric_face_flux",
               "Either 'volumetric_face_flux' or 'rhie_chow_user_object' must be provided.");
  const unsigned int compression_parameter_count = (_material_compression_flux ? 1 : 0) +
                                                   (_material_density ? 1 : 0) +
                                                   (_background_density ? 1 : 0);
  if (compression_parameter_count != 0 && compression_parameter_count != 3)
    paramError("material_compression_flux",
               "'material_compression_flux', 'material_density', and 'background_density' must "
               "be supplied together.");
}

Real
LinearFVLowMachMassAdvection::computeElemMatrixContribution()
{
  return _density_interp_coeffs.first * _face_volumetric_flux * _current_face_area;
}

Real
LinearFVLowMachMassAdvection::computeNeighborMatrixContribution()
{
  return _density_interp_coeffs.second * _face_volumetric_flux * _current_face_area;
}

Real
LinearFVLowMachMassAdvection::computeElemRightHandSideContribution()
{
  return -_face_explicit_correction * _current_face_area;
}

Real
LinearFVLowMachMassAdvection::computeNeighborRightHandSideContribution()
{
  return _face_explicit_correction * _current_face_area;
}

Real
LinearFVLowMachMassAdvection::computeBoundaryMatrixContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto & advective_bc = static_cast<const LinearFVAdvectionDiffusionBC &>(bc);
  const Real normal_factor = _current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0;
  return advective_bc.computeBoundaryValueMatrixContribution() * normal_factor *
         _face_volumetric_flux * _current_face_area;
}

Real
LinearFVLowMachMassAdvection::computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc)
{
  const auto & advective_bc = static_cast<const LinearFVAdvectionDiffusionBC &>(bc);
  const Real normal_factor = _current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0;
  return -advective_bc.computeBoundaryValueRHSContribution() * normal_factor *
         _face_volumetric_flux * _current_face_area;
}

void
LinearFVLowMachMassAdvection::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  _face_volumetric_flux = volumetricFaceFlux(*face_info);
  if (_current_face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    _density_interp_coeffs =
        interpCoeffs(Moose::FV::InterpMethod::Upwind, *face_info, true, _face_volumetric_flux);
    _face_explicit_correction = 0.0;
    if (_material_compression_flux)
    {
      const auto face_arg = functorFaceArg(*_material_compression_flux, *face_info);
      const Real compression_flux = (*_material_compression_flux)(face_arg, determineState());
      const Real material_density =
          (*_material_density)(functorFaceArg(*_material_density, *face_info), determineState());
      const Real background_density = (*_background_density)(
          functorFaceArg(*_background_density, *face_info), determineState());
      _face_explicit_correction = compression_flux * (material_density - background_density);
    }
  }
  else
  {
    _density_interp_coeffs = {0.0, 0.0};
    _face_explicit_correction = 0.0;
  }

  _assembled_explicit_correction[face_info->id()] = _face_explicit_correction;
}

Real
LinearFVLowMachMassAdvection::faceMassFluxDensity(const FaceInfo & face_info)
{
  const auto face_type = face_info.faceType(std::make_pair(_var_num, _sys_num));
  const auto state = determineState();
  const Real volumetric_flux = volumetricFaceFlux(face_info);

  if (face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    const auto correction = _assembled_explicit_correction.find(face_info.id());
    if (correction == _assembled_explicit_correction.end())
      mooseError(name(),
                 ": no assembled explicit correction is available for internal face ",
                 face_info.id(),
                 ". The mass-flux functor must be published after assembling the "
                 "temporary-density equation.");

    const Real donor_density = volumetric_flux >= 0.0
                                   ? _var.getElemValue(*face_info.elemInfo(), state)
                                   : _var.getElemValue(*face_info.neighborInfo(), state);
    return volumetric_flux * donor_density + correction->second;
  }

  if (face_type == FaceInfo::VarFaceNeighbors::ELEM ||
      face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR)
  {
    auto * const bc =
        dynamic_cast<LinearFVAdvectionDiffusionBC *>(_var.getBoundaryCondition(face_info));
    if (!bc)
    {
      if (std::abs(volumetric_flux) > TOLERANCE)
        mooseError(name(),
                   ": cannot publish nonzero boundary mass flux on face ",
                   face_info.id(),
                   " because the temporary-density variable has no compatible "
                   "LinearFVAdvectionDiffusionBC.");
      return 0.0;
    }

    bc->setupFaceData(&face_info, face_type);
    const auto * const elem_info = face_type == FaceInfo::VarFaceNeighbors::ELEM
                                       ? face_info.elemInfo()
                                       : face_info.neighborInfo();
    const Real cell_density = _var.getElemValue(*elem_info, state);
    const Real face_density = bc->computeBoundaryValueMatrixContribution() * cell_density +
                              bc->computeBoundaryValueRHSContribution();
    return volumetric_flux * face_density;
  }

  return 0.0;
}

Real
LinearFVLowMachMassAdvection::volumetricFaceFlux(const FaceInfo & face_info) const
{
  return _volumetric_face_flux
             ? (*_volumetric_face_flux)(functorFaceArg(*_volumetric_face_flux, face_info),
                                        determineState())
             : _volumetric_flux_provider->getVolumetricFaceFlux(face_info);
}
