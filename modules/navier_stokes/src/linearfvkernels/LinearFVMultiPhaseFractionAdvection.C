//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVMultiPhaseFractionAdvection.h"
#include "MooseLinearVariableFV.h"
#include "NSFVUtils.h"
#include "NS.h"
#include "LinearSystem.h"

registerMooseObject("NavierStokesApp", LinearFVMultiPhaseFractionAdvection);

InputParameters
LinearFVMultiPhaseFractionAdvection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription("Represents the matrix and right hand side contributions of an "
                             "advection term for a the mass-weighted phase fraction.");
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object",
      "The rhie-chow user-object which is used to determine the face velocity.");
  params.addRangeCheckedParam<Real>(
      "c_alpha", 0.0, "0.0<=c_alpha", "The compression velocity scaling constant.");
  params.addParam<MooseFunctorName>(NS::density, "The density.");
  params.addParam<bool>(
      "use_nonorthogonal_correction",
      false,
      "If the nonorthogonal correction should be used when computing the normal gradient.");
  params.addParam<bool>(
      "activate_mules",
      true,
      "Flag to aactivate CMULES limiting.");
  params.addParam<unsigned int>(
      "MULES_iterations",
      1,
      "Number of MULES iterations to perform.");
  params.addParam<bool>(
      "semi_implicit_mules",
      false,
      "Use semi-implicit MULES: two-pass approach with implicit upwind solve "
      "followed by explicitly limited anti-diffusive correction. Removes CFL "
      "restriction on the alpha equation.");
  params.addParam<bool>(
      "use_volumetric_flux",
      false,
      "Use volumetric face flux (phi) instead of mass flux (rho*phi) for the advection term. "
      "Required for purely volumetric alpha equations (no rho factor on time derivative). "
      "This matches the OpenFOAM interFoam formulation.");

  params += Moose::FV::advectedInterpolationParameter();

  MooseEnum limiterEnum(
      "min_mod vanLeer vanAlbada sou venkatakrishnan quick average upwind", "vanLeer");

  params.addParam<MooseEnum>("limiter_method",
                             limiterEnum,
                             "The limiter to use for the advected quantity. Options are "
                             "'min_mod', 'vanLeer', 'vanAlbada', 'sou', "
                             "'venkatakrishnan', 'quick', 'average', and "
                             "'upwind' with the default being 'upwind'.");
  return params;
}

LinearFVMultiPhaseFractionAdvection::LinearFVMultiPhaseFractionAdvection(
    const InputParameters & params)
  : LinearFVFluxKernel(params),
    _mass_flux_provider(const_cast<RhieChowMassFluxMultiPhase &>(getUserObject<RhieChowMassFluxMultiPhase>("rhie_chow_user_object"))),
    _dim(_subproblem.mesh().dimension()),
    _c_alpha(getParam<Real>("c_alpha")),
    _rho(params.isParamValid(NS::density) ? &(getFunctor<Real>(NS::density)) : nullptr),
    _use_nonorthogonal_correction(getParam<bool>("use_nonorthogonal_correction")),
    _use_mules(getParam<bool>("activate_mules")),
    _semi_implicit_mules(getParam<bool>("semi_implicit_mules")),
    _use_volumetric_flux(getParam<bool>("use_volumetric_flux")),
    _MULES_iterations(getParam<unsigned int>("MULES_iterations")),
    _advected_interp_coeffs(std::make_pair<Real, Real>(0, 0)),
    _total_adv_mass_face_flux(0.0),
    _limiter_method(getParam<MooseEnum>("limiter_method")),
    _compression_interp_coeffs(std::make_pair<Real, Real>(0, 0)),
    _total_comp_mass_face_flux(0.0)
{
  Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");

  if (_c_alpha > 1e-42)
  {
    if (!_rho && !_use_volumetric_flux)
      paramError(NS::density,
                 "The density must be provided when compression velocity is activated by setting "
                 "c_alpha>1e-42 (unless use_volumetric_flux=true).");

    // Gradients are needed for compression velocity
    _var.computeCellGradients();
  }

  // Semi-implicit MULES needs gradients for the post-solve correction step
  if (_semi_implicit_mules)
    _var.computeCellGradients();
}

Real
LinearFVMultiPhaseFractionAdvection::computeElemMatrixContribution()
{
  Real comp_mass_flux = 0.0;
  if (_c_alpha > 1e-42)
    comp_mass_flux += _compression_interp_coeffs.first * _total_comp_mass_face_flux;

  return (_advected_interp_coeffs.first * _total_adv_mass_face_flux + comp_mass_flux) *
         _current_face_area;
}

Real
LinearFVMultiPhaseFractionAdvection::computeNeighborMatrixContribution()
{
  Real comp_mass_flux = 0.0;
  if (_c_alpha > 1e-42)
    comp_mass_flux += _compression_interp_coeffs.second * _total_comp_mass_face_flux;

  return (_advected_interp_coeffs.second * _total_adv_mass_face_flux + comp_mass_flux) *
         _current_face_area;
}

Real
LinearFVMultiPhaseFractionAdvection::computeElemRightHandSideContribution()
{
  Real alpha_holo = this->getHighOrderFaceValue(_var) - this->getLowOrderFaceValue(_var);

  Real comp_mass_flux = 0.0;
  if(_c_alpha > 1e-42)
  {
    comp_mass_flux += _compression_interp_coeffs.first * _total_comp_mass_face_flux;
    if (_dim > 1 && _use_nonorthogonal_correction)
      comp_mass_flux += _compression_interp_coeffs.first *
                        computeCompressionVelocityMassFluxNonOrthogonalRHSContribution();
  }

  Real tol_mass_flux = _advected_interp_coeffs.first * _total_adv_mass_face_flux + comp_mass_flux;
  const auto rhs = _lambda_f * tol_mass_flux * alpha_holo * _current_face_area;

  return rhs;
}

Real
LinearFVMultiPhaseFractionAdvection::computeNeighborRightHandSideContribution()
{
  Real alpha_holo = this->getHighOrderFaceValue(_var) - this->getLowOrderFaceValue(_var);

  Real comp_mass_flux = 0.0;
  if(_c_alpha > 1e-42)
  {
    comp_mass_flux += _compression_interp_coeffs.second * _total_comp_mass_face_flux;
    if (_dim > 1 && _use_nonorthogonal_correction)
      comp_mass_flux += _compression_interp_coeffs.second *
                        computeCompressionVelocityMassFluxNonOrthogonalRHSContribution();
  }

  Real tol_mass_flux = -_advected_interp_coeffs.second * _total_adv_mass_face_flux - comp_mass_flux;
  const auto rhs = _lambda_f * tol_mass_flux * alpha_holo * _current_face_area;

  return rhs;
}

Real
LinearFVMultiPhaseFractionAdvection::computeBoundaryMatrixContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  const auto boundary_value_matrix_contrib = adv_bc->computeBoundaryValueMatrixContribution();

  // We support internal boundaries too so we have to make sure the normal points always outward
  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  // Adding compression velocity contribution
  Real compression_mass_flux = 0.0;
  if (_c_alpha > 1e-42)
    compression_mass_flux = computeCompressionVelocityMassFlux();

  return boundary_value_matrix_contrib * factor *
         (_total_adv_mass_face_flux + compression_mass_flux) * _current_face_area;
}

Real
LinearFVMultiPhaseFractionAdvection::computeBoundaryRHSContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  // We support internal boundaries too so we have to make sure the normal points always outward
  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0);

  const auto boundary_value_rhs_contrib = adv_bc->computeBoundaryValueRHSContribution();

  // Adding compression velocity contribution
  Real compression_mass_flux = 0.0;
  if (_c_alpha > 1e-42)
    compression_mass_flux = computeCompressionVelocityMassFlux();

  return -boundary_value_rhs_contrib * factor *
         (_total_adv_mass_face_flux + compression_mass_flux) * _current_face_area;
}

void
LinearFVMultiPhaseFractionAdvection::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  // Caching the flux on the face which will be reused in the advection term's matrix and RHS.
  // Volumetric mode uses phi (for purely volumetric alpha equations, matching OpenFOAM interFoam).
  // Mass mode uses rho*phi (for density-weighted alpha equations).
  if (_use_volumetric_flux)
    _total_adv_mass_face_flux = _mass_flux_provider.getVolumetricFaceFlux(*face_info);
  else
    _total_adv_mass_face_flux = _mass_flux_provider.getUnweightedMassFlux(*face_info);

  // Caching the interpolation coefficients so they will be reused for the matrix and right hand
  // side terms
  _advected_interp_coeffs =
      interpCoeffs(_advected_interp_method, *_current_face_info, true, _total_adv_mass_face_flux > 0);

  // Store low-order face
  _low_order_face = makeFace(
      *_current_face_info, limiterType(_advected_interp_method), _total_adv_mass_face_flux > 0);

  /// Get fluxes and interpolation cofficients for compression velocity
  if (_c_alpha > 1e-42)
  {
    _total_comp_mass_face_flux = computeCompressionVelocityMassFlux();
    _compression_interp_coeffs = 
      interpCoeffs(_advected_interp_method, *_current_face_info, true, _total_comp_mass_face_flux > 0);
  }

  // MULES
  if(_use_mules && !_semi_implicit_mules)
  {
    // Standard (explicit) MULES: compute lambda from old alpha field
    const auto total_adv_volume_flux =
        _mass_flux_provider.getVolumetricFaceFlux(*face_info) * _current_face_area * _dt / static_cast<Real>(_MULES_iterations);

    const bool donor_is_elem = _low_order_face.elem_is_upwind;

    auto donor =
        donor_is_elem ? _low_order_face.makeElem() : _low_order_face.makeNeighbor();
    auto acceptor =
        donor_is_elem ? _low_order_face.makeNeighbor() : _low_order_face.makeElem();

    auto donor_info = donor_is_elem ? _current_face_info->elemInfo()
                                    : _current_face_info->neighborInfo();
    auto acceptor_info = donor_is_elem ? _current_face_info->neighborInfo()
                                      : _current_face_info->elemInfo();

    if (!donor_info)
    {
      donor_info = acceptor_info;
      donor = acceptor;
    }

    if (!acceptor_info)
    {
      acceptor_info = donor_info;
      acceptor = donor;
    }

    const auto donnor_capacity =
        MetaPhysicL::raw_value(_var(donor, determineState())) * donor_info->volume();
    const auto acceptor_capacity =
        (1.0 - MetaPhysicL::raw_value(_var(acceptor, determineState()))) * acceptor_info->volume();

    _lambda_f = std::max(std::min(std::min(1.0, donnor_capacity / std::abs(total_adv_volume_flux)),
                                  acceptor_capacity / std::abs(total_adv_volume_flux)),
                        1e-10);
  }
  else if (_semi_implicit_mules)
  {
    // Semi-implicit MULES: assemble pure upwind system (lambda=0).
    // The correction is applied post-solve in applyMULESCorrection().
    _lambda_f = 0.0;
  }
  else
    _lambda_f = 1.0;

  // Store MULES-limited face alpha for consistent mass-momentum transport
  {
    const Real alpha_LO = this->getLowOrderFaceValue(_var);
    const Real alpha_HO = this->getHighOrderFaceValue(_var);
    const Real alpha_f_mules = alpha_LO + _lambda_f * (alpha_HO - alpha_LO);
    _mass_flux_provider.setConsistentFaceAlpha(*face_info, alpha_f_mules);
  }
}

Real
LinearFVMultiPhaseFractionAdvection::computeCompressionVelocityMassFluxNonOrthogonalRHSContribution()
{
  // Get the gradients from the adjacent cells
  const auto grad_elem = _var.gradSln(*_current_face_info->elemInfo());
  const auto & grad_neighbor = _var.gradSln(*_current_face_info->neighborInfo());

  // Interpolate the two gradients to the face
  const auto interp_coeffs =
      interpCoeffs(Moose::FV::InterpMethod::Average, *_current_face_info, true);

  const auto correction_vector =
      _current_face_info->normal() -
      1 / (_current_face_info->normal() * _current_face_info->eCN()) * _current_face_info->eCN();

  //  Compute compression velocity
  Real compression_velocity = computeCompressionVelocityMassFlux();

  return compression_velocity *
         (interp_coeffs.first * grad_elem + interp_coeffs.second * grad_neighbor) *
         correction_vector;
}

Real
LinearFVMultiPhaseFractionAdvection::computeCompressionVelocityMassFlux()
{

  // const auto alpha_f = this->getHighOrderFaceValue(_var);
  const auto alpha_f = this->getLowOrderFaceValue(_var);

  if (alpha_f <= 1e-12 || alpha_f >= 1.0 - 1e-12) // pure phase ⇒ no compression
    return 0.0;

  const auto grad = _var.gradSln(*_current_face_info->elemInfo());
  const auto grad_mag = grad.norm();

  if (grad_mag < 1e-14)
    return 0.0;

  const Real Un = std::fabs(_total_adv_mass_face_flux);

  const Real u_c = _c_alpha * Un;

  const auto compression_dir = (grad / grad_mag) * _current_face_info->normal();

  // In volumetric mode, compression flux is purely volumetric (no rho).
  // In mass mode, compression flux is density-weighted.
  if (_use_volumetric_flux)
  {
    return u_c * alpha_f * (1.0 - alpha_f) * compression_dir;
  }
  else
  {
    const auto rho = (*_rho)(_low_order_face, determineState());
    return rho * u_c * alpha_f * (1.0 - alpha_f) * compression_dir;
  }
}

Real
LinearFVMultiPhaseFractionAdvection::getLowOrderFaceValue(MooseLinearVariableFV<Real> & variable)
{
  return MetaPhysicL::raw_value(variable(_low_order_face, determineState()));
}

Real
LinearFVMultiPhaseFractionAdvection::getHighOrderFaceValue(MooseLinearVariableFV<Real> & variable)
{

  //---------------------------------------------------------------------------
  // 1. Donor / acceptor bookkeeping
  //---------------------------------------------------------------------------
  const bool donor_is_elem = _low_order_face.elem_is_upwind;

  auto donor    = donor_is_elem ? _low_order_face.makeElem()
                                : _low_order_face.makeNeighbor();
  auto acceptor = donor_is_elem ? _low_order_face.makeNeighbor() 
                                : _low_order_face.makeElem();

  const auto * donor_info    = donor_is_elem ? _current_face_info->elemInfo()
                                             : _current_face_info->neighborInfo();
  const auto * acceptor_info = donor_is_elem ? _current_face_info->neighborInfo()
                                             : _current_face_info->elemInfo();

  // Handle boundaries where one side is missing
  if (!donor_info)
  {
    donor_info = acceptor_info;
    donor = acceptor;
  }
  if (!acceptor_info)
  {
    acceptor_info = donor_info;
    acceptor = donor;
  }
    
  //---------------------------------------------------------------------------
  // 2. Cell-centred values and donor increment \Delta \phi_P (=\nabla \phi_P \cdot dP)
  //---------------------------------------------------------------------------
  const Real phi_P  = MetaPhysicL::raw_value(variable(donor, determineState()));
  const Real phi_N  = MetaPhysicL::raw_value(variable(acceptor, determineState()));
  const auto  gradP = variable.gradSln(*donor_info);               // \phi_P
  const Point face_c   = _current_face_info->faceCentroid();
  const Point donor_c  = donor_info->centroid();
  const auto  dP       = face_c - donor_c;                         // dP
  const Real delta_P   = gradP * dP;                               // \nabla \phi_P \cdot dP
  
  //---------------------------------------------------------------------------
  // 3. Slope ratio  r  and limiter \psi(r)
  //---------------------------------------------------------------------------

  constexpr Real tiny = 1.0e-14;

  const Real deltaPhi = phi_N - phi_P;
  const Real r        = delta_P / (deltaPhi + (deltaPhi >= 0 ? tiny : -tiny));

  Real psi = 1.0;   // default = second-order (\psi=1)

  // Map the MooseEnum value to the LimiterMethod enumerator
  LimiterMethod limiter_method = this->getLimiterMethod(_limiter_method);

  // Use the enumerator instead of string comparisons
  switch (limiter_method)
  {
    case MIN_MOD:
      psi = std::max(0.0, std::min(1.0, r));
      break;
    case VANLEER:
      psi = (r + std::fabs(r)) / (1.0 + std::fabs(r));
      break;
    case VANALBADA:
      psi = (r * r + r) / (r * r + 1.0);
      break;
    case QUICK:           // Koren QUICK
      psi = std::max(0.0,
                     std::min({ 2.0 / 3.0 * r + 1.0 / 6.0,  // bounded cubic
                                2.0 / 3.0,                  // upper plateaux
                                r }));                      // monotone
      break;
    case VENKATAKRISHNAN:
      psi = (r * r + 2.0 * r) / (r * r + r + 2.0);
      break;
    case AVERAGE:
      mooseError("`average` limiting not implemented. Please consider switching to average in the interpolation method.");
      break;
    case UPWIND:      // retain first-order
      psi = 0.0;
      break;
    default:
      psi = 1.0;  // 'sou' (second-order upwind) and anything unrecognised fall back to \psi = 1
      break;
  }

  //---------------------------------------------------------------------------
  // 4. High-order face value \phi_f  and storage for later access
  //---------------------------------------------------------------------------
  const Real phi_f = phi_P + psi * delta_P;  // Eq.  \phi_f = \phi_P + \psi(r)·\Delta \phi_P
  return phi_f;
}

void
LinearFVMultiPhaseFractionAdvection::applyMULESCorrection(LinearSystem & system)
{
  // Second pass of semi-implicit MULES.
  // The linear system was solved with lambda_f=0 (pure upwind), giving alpha_upwind.
  // Now compute bounded anti-diffusive correction using MULES limiters.
  //
  // All quantities are VOLUMETRIC (matching OpenFOAM interFoam):
  //   - phi = volumetric face flux [m^3/s]
  //   - capacity = alpha * V [m^3]
  //   - correction = dt/V * sum_f(lambda * phi * A * delta_alpha) [dimensionless]

  LinearImplicitSystem & li_system =
      libMesh::cast_ref<LinearImplicitSystem &>(system.system());
  NumericVector<Number> & solution = *(li_system.solution);

  // Update variable state and gradients to reflect alpha_upwind
  system.setSolution(*(li_system.current_local_solution));
  system.computeGradients();

  // Per-cell correction accumulator
  std::unordered_map<dof_id_type, Real> cell_correction;

  const auto & face_info_range = _subproblem.mesh().faceInfo();
  for (const auto * fi_ptr : face_info_range)
  {
    const auto & fi = *fi_ptr;

    if (!fi.neighborPtr())
      continue;

    const auto * elem_info = fi.elemInfo();
    const auto * neighbor_info = fi.neighborInfo();
    if (!elem_info || !neighbor_info)
      continue;

    // Always use VOLUMETRIC flux for the correction — this ensures dimensional
    // consistency: capacity [alpha*V] / flux [phi*A*dt] is dimensionless.
    const Real phi_vol = _mass_flux_provider.getVolumetricFaceFlux(fi);
    const Real face_area = fi.faceArea() * fi.faceCoord();

    // Upwind direction based on volumetric flux
    const bool elem_is_upwind = (phi_vol >= 0.0);

    Moose::FaceArg lo_face{&fi,
                           limiterType(_advected_interp_method),
                           elem_is_upwind,
                           false,
                           fi.elemPtr(),
                           nullptr};

    const Real alpha_LO = MetaPhysicL::raw_value(_var(lo_face, determineState()));

    // Donor/acceptor bookkeeping
    const bool donor_is_elem = elem_is_upwind;
    auto donor = donor_is_elem ? lo_face.makeElem() : lo_face.makeNeighbor();
    auto acceptor = donor_is_elem ? lo_face.makeNeighbor() : lo_face.makeElem();
    const auto * donor_info = donor_is_elem ? elem_info : neighbor_info;
    const auto * acceptor_info = donor_is_elem ? neighbor_info : elem_info;

    const Real phi_P = MetaPhysicL::raw_value(_var(donor, determineState()));
    const Real phi_N = MetaPhysicL::raw_value(_var(acceptor, determineState()));
    const auto gradP = _var.gradSln(*donor_info);
    const auto dP = fi.faceCentroid() - donor_info->centroid();
    const Real delta_P = gradP * dP;

    // Slope ratio and limiter
    constexpr Real tiny = 1.0e-14;
    const Real deltaPhi = phi_N - phi_P;
    const Real r = delta_P / (deltaPhi + (deltaPhi >= 0 ? tiny : -tiny));

    Real psi = 1.0;
    LimiterMethod lm = this->getLimiterMethod(_limiter_method);
    switch (lm)
    {
      case MIN_MOD:      psi = std::max(0.0, std::min(1.0, r)); break;
      case VANLEER:      psi = (r + std::fabs(r)) / (1.0 + std::fabs(r)); break;
      case VANALBADA:    psi = (r * r + r) / (r * r + 1.0); break;
      case QUICK:        psi = std::max(0.0, std::min({2.0/3.0*r + 1.0/6.0, 2.0/3.0, r})); break;
      case VENKATAKRISHNAN: psi = (r*r + 2.0*r) / (r*r + r + 2.0); break;
      case UPWIND:       psi = 0.0; break;
      default:           psi = 1.0; break;
    }

    const Real alpha_HO = phi_P + psi * delta_P;
    const Real alpha_diff = alpha_HO - alpha_LO;

    // Volumetric correction flux: phi_vol * A * alpha_diff [m^3/s]
    const Real phi_corr_vol = phi_vol * face_area * alpha_diff;

    // MULES limiter: capacity [m^3] / (flux * dt) [m^3] — dimensionless
    const Real donor_capacity = phi_P * donor_info->volume();
    const Real acceptor_capacity = (1.0 - phi_N) * acceptor_info->volume();
    const Real abs_corr_vol_dt = std::abs(phi_corr_vol * _dt);

    Real lambda_f = 1.0;
    if (abs_corr_vol_dt > 1e-42)
      lambda_f = std::max(std::min(std::min(1.0, donor_capacity / abs_corr_vol_dt),
                                    acceptor_capacity / abs_corr_vol_dt),
                          0.0);

    // Correction: delta_alpha = dt/V * lambda * phi_vol * A * alpha_diff
    const Real limited_corr = lambda_f * phi_corr_vol * _dt;

    const auto sys_num = _var.sys().number();
    const auto var_num = _var.number();
    const auto elem_dof = elem_info->dofIndices()[sys_num][var_num];
    const auto neighbor_dof = neighbor_info->dofIndices()[sys_num][var_num];

    cell_correction[elem_dof] -= limited_corr / elem_info->volume();
    cell_correction[neighbor_dof] += limited_corr / neighbor_info->volume();

    // Store MULES-limited face alpha for consistent mass-momentum transport
    const Real alpha_f_mules = alpha_LO + lambda_f * alpha_diff;
    _mass_flux_provider.setConsistentFaceAlpha(fi, alpha_f_mules);
  }

  // Apply corrections to the solution vector
  for (const auto & [dof_id, correction] : cell_correction)
    solution.add(dof_id, correction);

  solution.close();
  li_system.update();
  system.setSolution(*(li_system.current_local_solution));
}