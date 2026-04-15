//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVPLICAdvection.h"
#include "PLICReconstruction.h"
#include "PLICUtils3D.h"
#include "ElemInfo.h"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/elem.h"

registerMooseObject("NavierStokesApp", LinearFVPLICAdvection);

InputParameters
LinearFVPLICAdvection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Geometric VOF advection using PLIC interface reconstruction. "
      "Computes exact volume fluxes via face-plane clipping and Simpson's quadrature. "
      "Drop-in replacement for LinearFVMultiPhaseFractionAdvection.");
  params.addRequiredParam<UserObjectName>("plic_reconstruction",
                                          "The PLICReconstruction user object.");
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object",
      "The rhie-chow user-object which provides the volumetric face flux.");
  params.addParam<Real>(
      "alpha_tolerance", 1e-6, "Cells with alpha in (tol, 1-tol) use geometric flux.");
  return params;
}

LinearFVPLICAdvection::LinearFVPLICAdvection(const InputParameters & params)
  : LinearFVFluxKernel(params),
    _plic(getUserObject<PLICReconstruction>("plic_reconstruction")),
    _mass_flux_provider(getUserObject<RhieChowMassFluxMultiPhase>("rhie_chow_user_object")),
    _alpha_tol(getParam<Real>("alpha_tolerance")),
    _geometric_flux(0.0)
{
}

void
LinearFVPLICAdvection::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  _geometric_flux = 0.0;

  // Volumetric face flux: Q_f = u . n_face (m/s), positive = elem→neighbor
  const Real Q_f = _mass_flux_provider.getVolumetricFaceFlux(*face_info);

  if (std::abs(Q_f) < 1e-30)
    return;

  // Determine upwind cell
  const bool upwind_is_elem = (Q_f > 0.0);
  const ElemInfo * upwind_ei =
      upwind_is_elem ? face_info->elemInfo() : face_info->neighborInfo();

  if (!upwind_ei)
    return;

  // Read alpha from the current solution
  const auto & alpha_sys =
      libMesh::cast_ref<const libMesh::LinearImplicitSystem &>(_var.sys().system());
  const auto & solution = *alpha_sys.solution;
  const auto alpha_dof = upwind_ei->dofIndices()[_sys_num][_var_num];
  const Real alpha_up = solution(alpha_dof);

  // Full cell — all liquid flows through
  if (alpha_up >= 1.0 - _alpha_tol)
  {
    _geometric_flux = Q_f * _current_face_area * _dt;
    return;
  }

  // Empty cell — no liquid flows through
  if (alpha_up <= _alpha_tol)
    return;

  // Interfacial cell — need PLIC plane
  const auto upwind_elem_id = upwind_ei->elem()->id();
  if (!_plic.hasPlane(upwind_elem_id))
  {
    // No plane available (e.g. zero gradient) — fall back to upwind
    _geometric_flux = alpha_up * Q_f * _current_face_area * _dt;
    return;
  }

  const auto & plane = _plic.getPlane(upwind_elem_id);
  const auto & n_hat = plane.n_hat;
  const Real d = plane.d;

  // Extract face polygon vertices
  std::vector<Point> face_verts;
  getFaceVertices(*face_info, face_verts);

  // Interface-normal velocity component at the face:
  //   u_face = (Q_f / A_f) * n_face  (velocity vector at face)
  //   u_n = u_face . n_hat = (Q_f / A_f) * (n_face . n_hat)
  const Real u_n =
      (Q_f / face_info->faceArea()) * (face_info->normal() * n_hat);

  // Advect PLIC plane to three time levels
  const Real d_0 = d;
  const Real d_half = d - u_n * _dt * 0.5;
  const Real d_1 = d - u_n * _dt;

  // Simpson's rule time-integrated submerged area (m^2 * s)
  const Real time_integrated_area =
      NS::PLIC::timeIntegratedSubmergedArea(face_verts, n_hat, d_0, d_half, d_1, _dt);

  // Volume flux (m^3): Q_f (m/s) * time_integrated_area (m^2 * s)
  _geometric_flux = Q_f * time_integrated_area;
}

Real
LinearFVPLICAdvection::computeElemMatrixContribution()
{
  return 0.0;
}

Real
LinearFVPLICAdvection::computeNeighborMatrixContribution()
{
  return 0.0;
}

Real
LinearFVPLICAdvection::computeElemRightHandSideContribution()
{
  // Flux leaving elem reduces elem's alpha: RHS -= F_f
  return -_geometric_flux;
}

Real
LinearFVPLICAdvection::computeNeighborRightHandSideContribution()
{
  // Flux entering neighbor increases neighbor's alpha: RHS += F_f
  return _geometric_flux;
}

Real
LinearFVPLICAdvection::computeBoundaryMatrixContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "Boundary condition should be a LinearFVAdvectionDiffusionBC!");

  const auto boundary_value_matrix_contrib = adv_bc->computeBoundaryValueMatrixContribution();

  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  const Real Q_f = _mass_flux_provider.getVolumetricFaceFlux(*_current_face_info);

  return boundary_value_matrix_contrib * factor * Q_f * _current_face_area * _dt;
}

Real
LinearFVPLICAdvection::computeBoundaryRHSContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "Boundary condition should be a LinearFVAdvectionDiffusionBC!");

  const auto boundary_value_rhs_contrib = adv_bc->computeBoundaryValueRHSContribution();

  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  const Real Q_f = _mass_flux_provider.getVolumetricFaceFlux(*_current_face_info);

  return -boundary_value_rhs_contrib * factor * Q_f * _current_face_area * _dt;
}

void
LinearFVPLICAdvection::getFaceVertices(const FaceInfo & fi,
                                       std::vector<Point> & verts) const
{
  const auto & elem = fi.elem();
  const unsigned int side = fi.elemSideID();

  auto side_elem = elem.build_side_ptr(side);
  const unsigned int n_corners = side_elem->n_vertices();

  verts.resize(n_corners);
  for (unsigned int i = 0; i < n_corners; ++i)
    verts[i] = side_elem->point(i);
}
