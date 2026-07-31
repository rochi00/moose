//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LinearFVMaterialAdvection.h"

#include "LinearFVAdvectionDiffusionBC.h"
#include "GradientLimiterType.h"
#include "LowMachCUI.h"
#include "NSFVUtils.h"
#include "RhieChowMassFlux.h"

#include <cmath>

registerMooseObject("NavierStokesApp", LinearFVMaterialAdvection);

InputParameters
LinearFVMaterialAdvection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Implicitly discretizes div(u q) - q div(u) from one volumetric face flux so a "
      "constant advected field is preserved for non-divergence-free velocity.");
  params.addParam<MooseFunctorName>(
      "volumetric_face_flux",
      "",
      "Volumetric face-flux density oriented with the FaceInfo normal.");
  params.addParam<UserObjectName>(
      "rhie_chow_user_object", "", "Rhie-Chow user object providing the volumetric face flux.");
  params.addParam<MooseFunctorName>("u", "Advecting velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "Advecting velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "Advecting velocity in the z direction.");
  params += Moose::FV::advectedInterpolationParameter();
  return params;
}

LinearFVMaterialAdvection::LinearFVMaterialAdvection(const InputParameters & params)
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
    _u(isParamValid("u") && !getParam<MooseFunctorName>("u").empty() ? &getFunctor<Real>("u")
                                                                     : nullptr),
    _v(isParamValid("v") && !getParam<MooseFunctorName>("v").empty() ? &getFunctor<Real>("v")
                                                                     : nullptr),
    _w(isParamValid("w") && !getParam<MooseFunctorName>("w").empty() ? &getFunctor<Real>("w")
                                                                     : nullptr),
    _advected_interp_coeffs(0.0, 0.0),
    _face_volumetric_flux(0.0),
    _face_deferred_correction(0.0)
{
  const unsigned int number_of_flux_sources =
      (_volumetric_flux_provider ? 1 : 0) + (_volumetric_face_flux ? 1 : 0) + (_u ? 1 : 0);
  if (number_of_flux_sources != 1)
    paramError("volumetric_face_flux",
               "Exactly one of 'volumetric_face_flux', 'rhie_chow_user_object', or 'u' must be "
               "provided.");
  if ((_v || _w) && !_u)
    paramError("u", "The x velocity 'u' is required when 'v' or 'w' is supplied.");

  Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");
  if (_advected_interp_method != Moose::FV::InterpMethod::Upwind &&
      _advected_interp_method != Moose::FV::InterpMethod::Average &&
      _advected_interp_method != Moose::FV::InterpMethod::QUICK &&
      _advected_interp_method != Moose::FV::InterpMethod::Venkatakrishnan)
    paramError("advected_interp_method",
               "LinearFVMaterialAdvection supports 'upwind', 'average', legacy bounded cubic "
               "'quick', and geometry-aware 'venkatakrishnan' MUSCL interpolation.");
  if (_advected_interp_method == Moose::FV::InterpMethod::Venkatakrishnan)
    _var.computeCellGradients(Moose::FV::GradientLimiterType::Venkatakrishnan);
  else if (_advected_interp_method == Moose::FV::InterpMethod::QUICK)
    _var.computeCellGradients();
}

bool
LinearFVMaterialAdvection::usesDeferredCorrection() const
{
  return _advected_interp_method == Moose::FV::InterpMethod::QUICK ||
         _advected_interp_method == Moose::FV::InterpMethod::Venkatakrishnan;
}

void
LinearFVMaterialAdvection::addMatrixContribution()
{
  const Real integrated_flux = _face_volumetric_flux * _current_face_area;

  if (_current_face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    _dof_indices(0) = _current_face_info->elemInfo()->dofIndices()[_sys_num][_var_num];
    _dof_indices(1) = _current_face_info->neighborInfo()->dofIndices()[_sys_num][_var_num];

    const Real elem_coefficient = _advected_interp_coeffs.first;
    const Real neighbor_coefficient = _advected_interp_coeffs.second;

    // Element row: F q_f - F q_E.
    if (hasBlocks(_current_face_info->elemInfo()->subdomain_id()))
    {
      _matrix_contribution(0, 0) = (elem_coefficient - 1.0) * integrated_flux;
      _matrix_contribution(0, 1) = neighbor_coefficient * integrated_flux;
    }

    // Neighbor row: -F q_f + F q_N.
    if (hasBlocks(_current_face_info->neighborInfo()->subdomain_id()))
    {
      _matrix_contribution(1, 0) = -elem_coefficient * integrated_flux;
      _matrix_contribution(1, 1) = (1.0 - neighbor_coefficient) * integrated_flux;
    }

    for (auto & matrix : _matrices)
      matrix->add_matrix(_matrix_contribution, _dof_indices.get_values());
    return;
  }

  if (_current_face_type != FaceInfo::VarFaceNeighbors::ELEM &&
      _current_face_type != FaceInfo::VarFaceNeighbors::NEIGHBOR)
    return;

  auto * const bc = _var.getBoundaryCondition(*_current_face_info);
  const Real outward_flux =
      (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0) * integrated_flux;

  if (!bc)
  {
    if (std::abs(outward_flux) > TOLERANCE)
      mooseError(name(),
                 ": a nonzero material-advection flux crosses boundary face ",
                 _current_face_info->id(),
                 " without a LinearFVAdvectionDiffusionBC.");
    return;
  }

  const auto * const advective_bc = dynamic_cast<const LinearFVAdvectionDiffusionBC *>(bc);
  if (!advective_bc)
    mooseError(name(),
               ": boundary face ",
               _current_face_info->id(),
               " requires a LinearFVAdvectionDiffusionBC.");

  LinearFVFluxKernel::addMatrixContribution();
}

void
LinearFVMaterialAdvection::addRightHandSideContribution()
{
  if (_current_face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    if (!usesDeferredCorrection())
      return;

    _dof_indices(0) = _current_face_info->elemInfo()->dofIndices()[_sys_num][_var_num];
    _dof_indices(1) = _current_face_info->neighborInfo()->dofIndices()[_sys_num][_var_num];

    if (hasBlocks(_current_face_info->elemInfo()->subdomain_id()))
      _rhs_contribution(0) = -_face_deferred_correction * _current_face_area;
    if (hasBlocks(_current_face_info->neighborInfo()->subdomain_id()))
      _rhs_contribution(1) = _face_deferred_correction * _current_face_area;

    for (auto & vector : _vectors)
      vector->add_vector(_rhs_contribution.get_values().data(), _dof_indices.get_values());
    return;
  }

  if (_current_face_type != FaceInfo::VarFaceNeighbors::ELEM &&
      _current_face_type != FaceInfo::VarFaceNeighbors::NEIGHBOR)
    return;

  const Real integrated_flux = _face_volumetric_flux * _current_face_area;
  const Real outward_flux =
      (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0) * integrated_flux;
  auto * const bc = _var.getBoundaryCondition(*_current_face_info);

  if (!bc)
  {
    if (std::abs(outward_flux) > TOLERANCE)
      mooseError(name(),
                 ": a nonzero material-advection flux crosses boundary face ",
                 _current_face_info->id(),
                 " without a LinearFVAdvectionDiffusionBC.");
    return;
  }

  const auto * const advective_bc = dynamic_cast<const LinearFVAdvectionDiffusionBC *>(bc);
  if (!advective_bc)
    mooseError(name(),
               ": boundary face ",
               _current_face_info->id(),
               " requires a LinearFVAdvectionDiffusionBC.");

  LinearFVFluxKernel::addRightHandSideContribution();
}

Real
LinearFVMaterialAdvection::computeElemMatrixContribution()
{
  return (_advected_interp_coeffs.first - 1.0) * _face_volumetric_flux * _current_face_area;
}

Real
LinearFVMaterialAdvection::computeNeighborMatrixContribution()
{
  return _advected_interp_coeffs.second * _face_volumetric_flux * _current_face_area;
}

Real
LinearFVMaterialAdvection::computeElemRightHandSideContribution()
{
  return -_face_deferred_correction * _current_face_area;
}

Real
LinearFVMaterialAdvection::computeNeighborRightHandSideContribution()
{
  return _face_deferred_correction * _current_face_area;
}

Real
LinearFVMaterialAdvection::computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc)
{
  const auto & advective_bc = dynamic_cast<const LinearFVAdvectionDiffusionBC &>(bc);
  const Real outward_flux = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0) *
                            _face_volumetric_flux * _current_face_area;
  return (advective_bc.computeBoundaryValueMatrixContribution() - 1.0) * outward_flux;
}

Real
LinearFVMaterialAdvection::computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc)
{
  const auto & advective_bc = dynamic_cast<const LinearFVAdvectionDiffusionBC &>(bc);
  const Real outward_flux = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0) *
                            _face_volumetric_flux * _current_face_area;
  return -advective_bc.computeBoundaryValueRHSContribution() * outward_flux;
}

void
LinearFVMaterialAdvection::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  _face_volumetric_flux = volumetricFaceFlux(*face_info);
  if (_current_face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    if (usesDeferredCorrection())
    {
      _advected_interp_coeffs =
          interpCoeffs(Moose::FV::InterpMethod::Upwind, *face_info, true, _face_volumetric_flux);
      const bool upwind_is_elem = _face_volumetric_flux >= 0.0;
      const auto * const upwind_info =
          upwind_is_elem ? face_info->elemInfo() : face_info->neighborInfo();
      const Real donor_value = _var.getElemValue(*upwind_info, determineState());
      const Real high_order_value =
          _advected_interp_method == Moose::FV::InterpMethod::Venkatakrishnan
              ? computeMUSCLFaceValue(*face_info, upwind_is_elem)
              : computeCUIFaceValue(*face_info, upwind_is_elem);
      _face_deferred_correction = _face_volumetric_flux * (high_order_value - donor_value);
    }
    else
    {
      _advected_interp_coeffs =
          interpCoeffs(_advected_interp_method, *face_info, true, _face_volumetric_flux);
      _face_deferred_correction = 0.0;
    }
  }
  else
  {
    _advected_interp_coeffs = {0.0, 0.0};
    _face_deferred_correction = 0.0;
  }
}

Real
LinearFVMaterialAdvection::volumetricFaceFlux(const FaceInfo & face_info) const
{
  if (_volumetric_face_flux)
    return (*_volumetric_face_flux)(functorFaceArg(*_volumetric_face_flux, face_info),
                                    determineState());
  if (_volumetric_flux_provider)
    return _volumetric_flux_provider->getVolumetricFaceFlux(face_info);

  const auto state = determineState();
  RealVectorValue velocity;
  velocity(0) = (*_u)(functorFaceArg(*_u, face_info), state);
  if (_v)
    velocity(1) = (*_v)(functorFaceArg(*_v, face_info), state);
  if (_w)
    velocity(2) = (*_w)(functorFaceArg(*_w, face_info), state);
  return velocity * face_info.normal();
}

Real
LinearFVMaterialAdvection::computeCUIFaceValue(const FaceInfo & face_info,
                                               const bool upwind_is_elem) const
{
  const auto state = determineState();
  const auto * const upwind_info = upwind_is_elem ? face_info.elemInfo() : face_info.neighborInfo();
  const auto * const downwind_info =
      upwind_is_elem ? face_info.neighborInfo() : face_info.elemInfo();
  const Real upwind_value = _var.getElemValue(*upwind_info, state);
  const Real downwind_value = _var.getElemValue(*downwind_info, state);

  Point upwind_to_downwind = face_info.dCN();
  if (!upwind_is_elem)
    upwind_to_downwind *= -1.0;
  const Real far_upwind_value =
      downwind_value - 2.0 * (_var.gradSln(*upwind_info, state) * upwind_to_downwind);

  return LowMachCUI::faceValue(upwind_value, far_upwind_value, downwind_value);
}

Real
LinearFVMaterialAdvection::computeMUSCLFaceValue(const FaceInfo & face_info,
                                                 const bool upwind_is_elem) const
{
  const auto state = determineState();
  const auto * const upwind_info = upwind_is_elem ? face_info.elemInfo() : face_info.neighborInfo();
  const Real upwind_value = _var.getElemValue(*upwind_info, state);
  const auto upwind_gradient =
      _var.gradSln(*upwind_info, state, Moose::FV::GradientLimiterType::Venkatakrishnan);

  return upwind_value + upwind_gradient * (face_info.faceCentroid() - upwind_info->centroid());
}
