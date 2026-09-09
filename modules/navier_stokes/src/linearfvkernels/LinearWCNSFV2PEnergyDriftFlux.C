//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearWCNSFV2PEnergyDriftFlux.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", LinearWCNSFV2PEnergyDriftFlux);

InputParameters
LinearWCNSFV2PEnergyDriftFlux::validParams()
{
  auto params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Implements the enthalpy carried by the relative motion of the phases in the two-phase "
      "mixture model, div(beta_d beta_c / rho_m * (cp_d - cp_c) * T * u_slip), on the left hand "
      "side of the mixture energy equation.");
  params.addRequiredParam<MooseFunctorName>("u_slip", "The slip velocity in the x direction.");
  params.addParam<MooseFunctorName>("v_slip", "The slip velocity in the y direction.");
  params.addParam<MooseFunctorName>("w_slip", "The slip velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>("rho_d", "Dispersed phase density.");
  params.addRequiredParam<MooseFunctorName>("rho_c", "Continuous phase density.");
  params.addRequiredParam<MooseFunctorName>("cp_d", "Dispersed phase specific heat.");
  params.addRequiredParam<MooseFunctorName>("cp_c", "Continuous phase specific heat.");
  params.addParam<MooseFunctorName>("fd", 0.0, "Fraction dispersed phase.");
  params.renameParam("fd", "fraction_dispersed", "");

  // This term advects the temperature, so it needs a boundary value to act on and must only run
  // where the variable carries a boundary condition. 'force_boundary_execution' is deliberately
  // left at its default of false: forcing it hands this kernel a null boundary condition on any
  // boundary the temperature is not constrained on.
  params.addParam<std::vector<BoundaryName>>(
      "slip_boundaries",
      {},
      "Boundaries across which the dispersed phase may travel, normally the inlets and the "
      "outlets. The enthalpy carried by the relative motion contributes on these boundaries and "
      "on every internal face, but not on impermeable boundaries such as walls, where the phase "
      "cannot cross even though the temperature may carry a boundary condition.");

  MooseEnum coeff_interp_method("average harmonic", "average");
  params.addParam<MooseEnum>(
      "coeff_interp_method",
      coeff_interp_method,
      "Switch that can select the face interpolation method for the enthalpy flux coefficient.");

  params += Moose::FV::advectedInterpolationParameter();
  return params;
}

LinearWCNSFV2PEnergyDriftFlux::LinearWCNSFV2PEnergyDriftFlux(const InputParameters & params)
  : LinearFVFluxKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _rho_d(getFunctor<Real>("rho_d")),
    _rho_c(getFunctor<Real>("rho_c")),
    _cp_d(getFunctor<Real>("cp_d")),
    _cp_c(getFunctor<Real>("cp_c")),
    _f_d(getFunctor<Real>("fd")),
    _u_slip(getFunctor<Real>("u_slip")),
    _v_slip(isParamValid("v_slip") ? &getFunctor<Real>("v_slip") : nullptr),
    _w_slip(isParamValid("w_slip") ? &getFunctor<Real>("w_slip") : nullptr),
    _coeff_interp_method(
        Moose::FV::selectInterpolationMethod(getParam<MooseEnum>("coeff_interp_method"))),
    _advected_interp_coeffs(std::make_pair<Real, Real>(0, 0)),
    _face_flux(0.0),
    _boundary_normal_factor(1.0)
{
  Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");

  for (const auto & bname : getParam<std::vector<BoundaryName>>("slip_boundaries"))
    _slip_boundaries.insert(_mesh.getBoundaryID(bname));

  if (_dim >= 2 && !_v_slip)
    mooseError("In two or more dimensions, the v_slip velocity must be supplied using the 'v_slip' "
               "parameter");
  if (_dim >= 3 && !_w_slip)
    mooseError(
        "In three dimensions, the w_slip velocity must be supplied using the 'w_slip' parameter");
}

Real
LinearWCNSFV2PEnergyDriftFlux::computeElemMatrixContribution()
{
  return _advected_interp_coeffs.first * _face_flux * _current_face_area;
}

Real
LinearWCNSFV2PEnergyDriftFlux::computeNeighborMatrixContribution()
{
  return _advected_interp_coeffs.second * _face_flux * _current_face_area;
}

Real
LinearWCNSFV2PEnergyDriftFlux::computeElemRightHandSideContribution()
{
  // The term is linear in the advected temperature, so it is carried entirely by the matrix
  return 0.0;
}

Real
LinearWCNSFV2PEnergyDriftFlux::computeNeighborRightHandSideContribution()
{
  return 0.0;
}

Real
LinearWCNSFV2PEnergyDriftFlux::computeBoundaryMatrixContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  const auto boundary_value_matrix_contrib = adv_bc->computeBoundaryValueMatrixContribution();

  return boundary_value_matrix_contrib * _boundary_normal_factor * _face_flux * _current_face_area;
}

Real
LinearWCNSFV2PEnergyDriftFlux::computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  const auto boundary_value_rhs_contrib = adv_bc->computeBoundaryValueRHSContribution();

  return -boundary_value_rhs_contrib * _boundary_normal_factor * _face_flux * _current_face_area;
}

void
LinearWCNSFV2PEnergyDriftFlux::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  // Multiplier that ensures the normal of the boundary always points outwards, even in cases
  // when the boundary is within the mesh
  _boundary_normal_factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  const auto & normal = _current_face_info->normal();
  const auto state = determineState();
  const bool on_boundary = Moose::FV::onBoundary(*this, *_current_face_info);

  const auto face_arg =
      on_boundary ? singleSidedFaceArg(_current_face_info) : makeCDFace(*_current_face_info);

  RealVectorValue u_slip_vel_vec;
  if (_dim == 1)
    u_slip_vel_vec = RealVectorValue(_u_slip(face_arg, state), 0.0, 0.0);
  else if (_dim == 2)
    u_slip_vel_vec = RealVectorValue(_u_slip(face_arg, state), (*_v_slip)(face_arg, state), 0.0);
  else
    u_slip_vel_vec = RealVectorValue(
        _u_slip(face_arg, state), (*_v_slip)(face_arg, state), (*_w_slip)(face_arg, state));

  Real face_coefficient;
  if (on_boundary)
    face_coefficient = enthalpyFluxCoefficient(face_arg, state);
  else
  {
    const auto elem_arg = makeElemArg(_current_face_info->elemPtr());
    const auto neigh_arg = makeElemArg(_current_face_info->neighborPtr());

    Moose::FV::interpolate(_coeff_interp_method,
                           face_coefficient,
                           enthalpyFluxCoefficient(elem_arg, state),
                           enthalpyFluxCoefficient(neigh_arg, state),
                           *_current_face_info,
                           true);
  }

  // Positive on the left hand side: a LinearFVFluxKernel assembles its face flux onto the left
  // hand side, and the phase-summed advection leaves +div(sum_k a_k rho_k h_k u_Mk) there.
  //
  // Internal faces always receive the flux. A boundary face receives it only if the boundary was
  // declared permeable, since the phases cannot separate across an impermeable wall and carrying
  // enthalpy through one would be a leak.
  const bool slip_allowed_here =
      _current_face_info->neighborPtr() ||
      (!_current_face_info->boundaryIDs().empty() &&
       _slip_boundaries.count(*_current_face_info->boundaryIDs().begin()));

  _face_flux = slip_allowed_here ? face_coefficient * (normal * u_slip_vel_vec) : 0.0;

  // Only internal faces need advected interpolation coefficients; boundary faces are handled
  // through the linear FV boundary conditions
  if (_current_face_type == FaceInfo::VarFaceNeighbors::BOTH)
    _advected_interp_coeffs =
        interpCoeffs(_advected_interp_method, *_current_face_info, true, _face_flux);
}
