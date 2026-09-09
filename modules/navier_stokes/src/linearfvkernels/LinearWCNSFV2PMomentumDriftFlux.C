//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearWCNSFV2PMomentumDriftFlux.h"
#include "NS.h"
#include "RhieChowMassFlux.h"
#include "LinearFVBoundaryCondition.h"
#include "LinearFVAdvectionDiffusionBC.h"

#include <limits>

registerMooseObject("NavierStokesApp", LinearWCNSFV2PMomentumDriftFlux);
// Renamed for consistency with the other linear finite volume Navier Stokes objects,
// which all carry the WCNSLinearFV prefix
registerMooseObjectRenamed("NavierStokesApp",
                           LinearWCNSFV2PMomentumDriftFlux,
                           "08/18/2027 00:00",
                           LinearWCNSFV2PMomentumDriftFlux);

InputParameters
LinearWCNSFV2PMomentumDriftFlux ::validParams()
{
  auto params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Implements the diffusion (drift) stress of the two-phase mixture model, "
      "div(beta_d beta_c / rho_m * u_slip (x) u_slip), on the left hand side of the mixture "
      "momentum equation.");
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object",
      "The rhie-chow user-object which is used to determine the face velocity.");
  params.addRequiredParam<MooseFunctorName>("u_slip", "The slip velocity in the x direction.");
  params.addParam<MooseFunctorName>("v_slip", "The slip velocity in the y direction.");
  params.addParam<MooseFunctorName>("w_slip", "The slip velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>("rho_d", "Dispersed phase density.");
  params.addRequiredParam<MooseFunctorName>("rho_c", "Continuous phase density.");
  params.addParam<MooseFunctorName>("fd", 0.0, "Fraction dispersed phase.");
  params.renameParam("fd", "fraction_dispersed", "");

  params.addParam<bool>(
      "force_boundary_execution", true, "This kernel should execute on boundaries by default");
  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>(
      "momentum_component",
      momentum_component,
      "The component of the momentum equation that this kernel applies to.");

  MooseEnum coeff_interp_method("average harmonic", "harmonic");
  params.addParam<MooseEnum>("density_interp_method",
                             coeff_interp_method,
                             "Switch that can select face interpolation method for the density.");

  return params;
}

LinearWCNSFV2PMomentumDriftFlux ::LinearWCNSFV2PMomentumDriftFlux(const InputParameters & params)
  : LinearFVFluxKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _mass_flux_provider(getUserObject<RhieChowMassFlux>("rhie_chow_user_object")),
    _rho_d(getFunctor<Real>("rho_d")),
    _rho_c(getFunctor<Real>("rho_c")),
    _f_d(getFunctor<Real>("fd")),
    _u_slip(getFunctor<Real>("u_slip")),
    _v_slip(isParamValid("v_slip") ? &getFunctor<Real>("v_slip") : nullptr),
    _w_slip(isParamValid("w_slip") ? &getFunctor<Real>("w_slip") : nullptr),
    _index(getParam<MooseEnum>("momentum_component")),
    _density_interp_method(
        Moose::FV::selectInterpolationMethod(getParam<MooseEnum>("density_interp_method"))),
    _face_flux(0.0),
    _slip_mass_flux(0.0),
    _gamma(0.0),
    _boundary_normal_factor(1.0)
{
  if (_dim >= 2 && !_v_slip)
    mooseError("In two or more dimensions, the v_slip velocity must be supplied using the 'v_slip' "
               "parameter");
  if (_dim >= 3 && !_w_slip)
    mooseError(
        "In three dimensions, the w_slip velocity must be supplied using the 'w_slip' parameter");

  // Note that since this is used in a segregated solver at this time, we don't need to declare
  // that this kernel may depend on a phase fraction variable
}

void
LinearWCNSFV2PMomentumDriftFlux::computeFlux()
{
  const auto & normal = _current_face_info->normal();
  const auto state = determineState();
  const bool on_boundary = Moose::FV::onBoundary(*this, *_current_face_info);

  Moose::FaceArg face_arg;
  if (on_boundary)
    face_arg = singleSidedFaceArg(_current_face_info);
  else
    face_arg = makeCDFace(*_current_face_info);

  RealVectorValue u_slip_vel_vec;
  if (_dim == 1)
    u_slip_vel_vec = RealVectorValue(_u_slip(face_arg, state), 0.0, 0.0);
  else if (_dim == 2)
    u_slip_vel_vec = RealVectorValue(_u_slip(face_arg, state), (*_v_slip)(face_arg, state), 0.0);
  else
    u_slip_vel_vec = RealVectorValue(
        _u_slip(face_arg, state), (*_v_slip)(face_arg, state), (*_w_slip)(face_arg, state));

  const auto uslipdotn = normal * u_slip_vel_vec;

  // The exact diffusion stress coefficient, beta_d beta_c / rho_m, evaluated on the face
  Real face_coefficient;
  if (on_boundary)
    face_coefficient = diffusionStressCoefficient(face_arg, state);
  else
  {
    const auto elem_arg = makeElemArg(_current_face_info->elemPtr());
    const auto neigh_arg = makeElemArg(_current_face_info->neighborPtr());

    const auto elem_coefficient = diffusionStressCoefficient(elem_arg, state);
    const auto neighbor_coefficient = diffusionStressCoefficient(neigh_arg, state);

    // beta_d beta_c / rho_m vanishes wherever either phase is absent, at a phase fraction of zero
    // and again at one, and the harmonic mean is not defined there. Fall back to the arithmetic
    // average on those faces. The coefficient weights an advective flux rather than acting as a
    // diffusivity, so the arithmetic average is the natural choice for it in any case; the
    // harmonic option is retained for continuity with the parameter's previous meaning.
    const auto interp_method = (elem_coefficient > 0.0 && neighbor_coefficient > 0.0)
                                   ? _density_interp_method
                                   : Moose::FV::InterpMethod::Average;

    Moose::FV::interpolate(interp_method,
                           face_coefficient,
                           elem_coefficient,
                           neighbor_coefficient,
                           *_current_face_info,
                           true);
  }

  // The term is written as a flux carried by the slip velocity, so that the flux scale can be
  // reused below as the scale of the implicit surrogate.
  //
  // Sign. A LinearFVFluxKernel assembles its face flux onto the left hand side, so the flux set
  // here is the left hand side form of the term. Summing the phase momentum equations puts
  // +div(sum_k a_k rho_k u_Mk u_Mk) on the left hand side, hence the positive sign. This differs
  // from the nonlinear WCNSFV2PMomentumDriftFlux, which carries the opposite sign and is left
  // unchanged; the two discretizations therefore disagree on this term by construction.
  _slip_mass_flux = face_coefficient * uslipdotn;
  _face_flux = _slip_mass_flux * u_slip_vel_vec(_index);
}

Real
LinearWCNSFV2PMomentumDriftFlux::computeElemMatrixContribution()
{
  return std::max(_gamma, 0.0) * _current_face_area;
}

Real
LinearWCNSFV2PMomentumDriftFlux::computeNeighborMatrixContribution()
{
  return -std::max(-_gamma, 0.0) * _current_face_area;
}

Real
LinearWCNSFV2PMomentumDriftFlux::deferredCorrection() const
{
  // Upwind value of the previous iterate with respect to the slip mass flux
  const auto old_state = Moose::previousNonlinearState();
  const auto u_upwind_old =
      _gamma > 0 ? _var(makeElemArg(_current_face_info->elemPtr()), old_state).value()
                 : _var(makeElemArg(_current_face_info->neighborPtr()), old_state).value();

  return (_gamma * u_upwind_old - _face_flux) * _current_face_area;
}

Real
LinearWCNSFV2PMomentumDriftFlux::computeElemRightHandSideContribution()
{
  return deferredCorrection();
}

Real
LinearWCNSFV2PMomentumDriftFlux::computeNeighborRightHandSideContribution()
{
  // The right hand side contributions are not negated by the assembly routine, unlike the matrix
  // ones, so the neighbour row has to be given the opposite sign here
  return -deferredCorrection();
}

Real
LinearWCNSFV2PMomentumDriftFlux::computeBoundaryRHSContribution(
    const LinearFVBoundaryCondition & /*bc*/)
{
  // Lagging the whole term for now
  // TODO: make sure this only gets called once, and not once per BC
  return -_boundary_normal_factor * _face_flux * _current_face_area;
}

void
LinearWCNSFV2PMomentumDriftFlux::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  // Multiplier that ensures the normal of the boundary always points outwards, even in cases
  // when the boundary is within the mesh
  _boundary_normal_factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  // Caching the flux on the face which will be reused in the matrix and right hand side
  // contributions
  computeFlux();

  // Coefficient of the implicit surrogate. This term has no genuine linear dependence on the
  // velocity component being solved for: the slip velocity is driven by the pressure gradient and
  // by gravity through its algebraic closure, not by the local cell value. The matrix contribution
  // is therefore a deferred correction rather than a linearization, and the deferredCorrection()
  // added to the right hand side cancels it exactly at convergence. The converged solution does
  // not depend on _gamma, which is chosen for convergence alone.
  //
  // The natural choice is the Picard ratio flux / u_old, but that grows without bound as u_old
  // approaches zero and carries no controlled sign. We cap its magnitude with the slip mass flux,
  // which is the scale of the term itself, and take the sign of the slip mass flux. Combined with
  // the upwind structure of the matrix contributions this puts a non-negative coefficient on both
  // diagonals and a non-positive one on both off-diagonals, on every face, with no special case
  // near stagnation.
  const auto u_old = _var(makeCDFace(*_current_face_info), Moose::previousNonlinearState()).value();
  const auto ratio = (std::abs(u_old) > libMesh::TOLERANCE)
                         ? std::abs(_face_flux / u_old)
                         : std::numeric_limits<Real>::max();
  _gamma = std::copysign(std::min(ratio, std::abs(_slip_mass_flux)), _slip_mass_flux);
}
