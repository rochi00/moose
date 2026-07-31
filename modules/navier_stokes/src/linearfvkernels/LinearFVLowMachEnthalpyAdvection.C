//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LinearFVLowMachEnthalpyAdvection.h"

#include "LowMachCUI.h"
#include "NSFVUtils.h"
#include "RhieChowMassFlux.h"

registerMooseObject("NavierStokesApp", LinearFVLowMachEnthalpyAdvection);

InputParameters
LinearFVLowMachEnthalpyAdvection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Advects fixed-outer-iterate low-Mach specific enthalpy with bounded CUI reconstruction and "
      "the common mass-equation face flux.");
  params.addRequiredParam<MooseFunctorName>(
      "specific_enthalpy", "Specific enthalpy fixed at the current outer iterate.");
  params.addRequiredParam<MooseFunctorName>(
      "dh_dT", "Current enthalpy derivative used to reconstruct the upwind enthalpy gradient.");
  params.addParam<MooseFunctorName>(
      "inflow_specific_enthalpy",
      "Specific enthalpy imposed where a boundary face carries mass into the domain.");
  params.addParam<UserObjectName>(
      "rhie_chow_user_object", "", "Rhie-Chow user object providing the common face mass flux.");
  params.addParam<MooseFunctorName>(
      "mass_flux_functor",
      "",
      "Optional common face mass-flux functor. This overrides the Rhie-Chow user object.");
  params.addParam<bool>(
      "mass_flux_is_integrated", false, "Whether mass_flux_functor already includes face area.");
  return params;
}

LinearFVLowMachEnthalpyAdvection::LinearFVLowMachEnthalpyAdvection(const InputParameters & params)
  : LinearFVFluxKernel(params),
    _specific_enthalpy(getFunctor<Real>("specific_enthalpy")),
    _dh_dT(getFunctor<Real>("dh_dT")),
    _inflow_specific_enthalpy(isParamValid("inflow_specific_enthalpy")
                                  ? &getFunctor<Real>("inflow_specific_enthalpy")
                                  : nullptr),
    _mass_flux_provider(isParamValid("rhie_chow_user_object") &&
                                !getParam<UserObjectName>("rhie_chow_user_object").empty()
                            ? &getUserObject<RhieChowMassFlux>("rhie_chow_user_object")
                            : nullptr),
    _mass_flux_functor(isParamValid("mass_flux_functor") &&
                               !getParam<MooseFunctorName>("mass_flux_functor").empty()
                           ? &getFunctor<Real>("mass_flux_functor")
                           : nullptr),
    _mass_flux_functor_name(getParam<MooseFunctorName>("mass_flux_functor")),
    _mass_flux_is_integrated(getParam<bool>("mass_flux_is_integrated")),
    _face_contribution(0.0)
{
  if (!_mass_flux_provider && !_mass_flux_functor)
    paramError("mass_flux_functor",
               "Either 'mass_flux_functor' or 'rhie_chow_user_object' must be provided.");

  _var.computeCellGradients();
}

Real
LinearFVLowMachEnthalpyAdvection::computeElemMatrixContribution()
{
  return 0.0;
}

Real
LinearFVLowMachEnthalpyAdvection::computeNeighborMatrixContribution()
{
  return 0.0;
}

Real
LinearFVLowMachEnthalpyAdvection::computeElemRightHandSideContribution()
{
  return -_face_contribution;
}

Real
LinearFVLowMachEnthalpyAdvection::computeNeighborRightHandSideContribution()
{
  return _face_contribution;
}

Real
LinearFVLowMachEnthalpyAdvection::computeBoundaryMatrixContribution(
    const LinearFVBoundaryCondition & /*bc*/)
{
  return 0.0;
}

Real
LinearFVLowMachEnthalpyAdvection::computeBoundaryRHSContribution(
    const LinearFVBoundaryCondition & /*bc*/)
{
  const Real normal_factor = _current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0;
  return -normal_factor * _face_contribution;
}

void
LinearFVLowMachEnthalpyAdvection::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  if (_use_fixed_point_flux)
  {
    const auto fixed_flux = _fixed_point_face_contribution.find(face_info->id());
    if (fixed_flux == _fixed_point_face_contribution.end())
      mooseError(
          name(), ": no fixed-point enthalpy flux was captured for face ", face_info->id(), ".");
    _face_contribution = fixed_flux->second;
  }
  else
  {
    const auto state = determineState();
    const bool internal_face = _current_face_type == FaceInfo::VarFaceNeighbors::BOTH;
    const Moose::FaceArg mass_flux_arg =
        internal_face ? makeCDFace(*face_info) : singleSidedFaceArg(face_info);

    Real mass_flux = _mass_flux_functor ? (*_mass_flux_functor)(mass_flux_arg, state)
                                        : _mass_flux_provider->getMassFlux(*face_info);
    if (_mass_flux_functor && _mass_flux_is_integrated)
      mass_flux = _current_face_area > 0.0 ? mass_flux / _current_face_area : 0.0;

    Real face_enthalpy;
    if (internal_face)
      face_enthalpy = computeCUIFaceEnthalpy(*face_info, mass_flux >= 0.0);
    else
    {
      const Real outward_normal =
          _current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0;
      if (outward_normal * mass_flux < 0.0)
      {
        if (!_inflow_specific_enthalpy)
          mooseError(name(),
                     ": boundary face ",
                     face_info->id(),
                     " carries mass into the domain, but no 'inflow_specific_enthalpy' was "
                     "provided.");
        face_enthalpy = (*_inflow_specific_enthalpy)(
            functorFaceArg(*_inflow_specific_enthalpy, *face_info), state);
      }
      else
        face_enthalpy = _specific_enthalpy(functorFaceArg(_specific_enthalpy, *face_info), state);
    }

    _face_contribution = mass_flux * face_enthalpy * _current_face_area;
    if (_capture_fixed_point_flux)
      _fixed_point_face_contribution[face_info->id()] = _face_contribution;
  }

  _assembled_face_contribution[face_info->id()] = _face_contribution;
}

void
LinearFVLowMachEnthalpyAdvection::beginFixedPointFluxCapture()
{
  _fixed_point_face_contribution.clear();
  _capture_fixed_point_flux = true;
  _use_fixed_point_flux = false;
}

void
LinearFVLowMachEnthalpyAdvection::finishFixedPointFluxCapture()
{
  _capture_fixed_point_flux = false;
  _use_fixed_point_flux = true;
}

Real
LinearFVLowMachEnthalpyAdvection::assembledFaceEnthalpyFlux(const FaceInfo & face_info) const
{
  const auto flux = _assembled_face_contribution.find(face_info.id());
  if (flux == _assembled_face_contribution.end())
    mooseError(name(),
               ": no assembled enthalpy flux is available for face ",
               face_info.id(),
               ". The low-Mach divergence source must be constructed after the final enthalpy "
               "equation assembly.");
  return flux->second;
}

Real
LinearFVLowMachEnthalpyAdvection::computeCUIFaceEnthalpy(const FaceInfo & face_info,
                                                         const bool upwind_is_elem) const
{
  const auto state = determineState();
  const auto * const upwind_info = upwind_is_elem ? face_info.elemInfo() : face_info.neighborInfo();
  const auto * const downwind_info =
      upwind_is_elem ? face_info.neighborInfo() : face_info.elemInfo();
  const auto upwind_arg = makeElemArg(upwind_info->elem());
  const auto downwind_arg = makeElemArg(downwind_info->elem());
  const Real upwind_enthalpy = _specific_enthalpy(upwind_arg, state);
  const Real downwind_enthalpy = _specific_enthalpy(downwind_arg, state);

  Point upwind_to_downwind = face_info.dCN();
  if (!upwind_is_elem)
    upwind_to_downwind *= -1.0;

  Real far_upwind_enthalpy;
  bool reconstructed_from_inflow = false;
  if (_inflow_specific_enthalpy)
    try
    {
      const unsigned int upwind_side =
          upwind_is_elem ? face_info.elemSideID() : face_info.neighborSideID();
      const unsigned int far_upwind_side = upwind_info->elem()->opposite_side(upwind_side);
      if (!upwind_info->elem()->neighbor_ptr(far_upwind_side))
      {
        const auto * const boundary_face = _mesh.faceInfo(upwind_info->elem(), far_upwind_side);
        const auto boundary_face_type = boundary_face->faceType(std::make_pair(_var_num, _sys_num));
        const Real outward_normal = boundary_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0
                                    : boundary_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR
                                        ? -1.0
                                        : 0.0;
        const Moose::FaceArg boundary_flux_arg = singleSidedFaceArg(boundary_face);
        const Real boundary_mass_flux = _mass_flux_functor
                                            ? (*_mass_flux_functor)(boundary_flux_arg, state)
                                            : _mass_flux_provider->getMassFlux(*boundary_face);
        const Real cell_spacing = upwind_to_downwind.norm();
        if (outward_normal * boundary_mass_flux < 0.0 && cell_spacing > TOLERANCE)
        {
          const Real boundary_distance =
              -((boundary_face->faceCentroid() - upwind_info->centroid()) *
                (upwind_to_downwind / cell_spacing));
          if (boundary_distance > TOLERANCE)
          {
            const Real boundary_enthalpy = (*_inflow_specific_enthalpy)(
                functorFaceArg(*_inflow_specific_enthalpy, *boundary_face), state);
            far_upwind_enthalpy = upwind_enthalpy + cell_spacing / boundary_distance *
                                                        (boundary_enthalpy - upwind_enthalpy);
            reconstructed_from_inflow = true;
          }
        }
      }
    }
    catch (libMesh::LogicError &)
    {
      // Non-tensor cells have no unique opposite side; use the multidimensional gradient below.
    }

  if (!reconstructed_from_inflow)
    far_upwind_enthalpy =
        downwind_enthalpy -
        2.0 * _dh_dT(upwind_arg, state) * (_var.gradSln(*upwind_info, state) * upwind_to_downwind);

  return LowMachCUI::faceValue(upwind_enthalpy, far_upwind_enthalpy, downwind_enthalpy);
}
