//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVFluxKernel.h"
#include "LinearFVAdvectionDiffusionBC.h"
#include "MathFVUtils.h"

#include <algorithm>

/**
 * Enthalpy carried by the relative motion of the phases, the energy counterpart of the diffusion
 * (drift) stress in the mixture momentum equation.
 *
 * The exact advective term of the mixture energy equation is the sum over the phases of each phase
 * carrying its own enthalpy at its own velocity,
 *
 * \f[
 *   \nabla \cdot \sum_k \alpha_k \rho_k h_k u_k
 * \f]
 *
 * Substituting \f$ u_k = u_m + u_{Mk} \f$ splits this into the mixture term that
 * LinearFVEnergyAdvection already assembles plus a remainder,
 *
 * \f[
 *   \nabla \cdot \left( \rho_m h_m u_m \right)
 *     + \nabla \cdot \sum_k \alpha_k \rho_k h_k u_{Mk}
 * \f]
 *
 * For a single dispersed phase, using \f$ \beta_c u_{Mc} = -\beta_d u_{Md} \f$ and
 * \f$ u_{Md} = (\beta_c / \rho_m) u_s \f$, the remainder collapses onto the slip velocity with the
 * same coefficient that carries the diffusion stress,
 *
 * \f[
 *   \sum_k \alpha_k \rho_k h_k u_{Mk}
 *     = \frac{\beta_d \beta_c}{\rho_m} \left( h_d - h_c \right) u_s
 * \f]
 *
 * With both phases at the mixture temperature, \f$ h_k = c_{p,k} T \f$, this kernel therefore
 * assembles \f$ \nabla \cdot \left[ \frac{\beta_d \beta_c}{\rho_m} (c_{p,d} - c_{p,c}) T u_s
 * \right] \f$ on the left hand side. It vanishes identically when the two specific heats are
 * equal.
 *
 * This is the term written \f$ q_p = \sum_k \alpha_k \rho_k e_k u_k \f$ in the INL reference
 * report, and the phase-summed advection of ANSYS Fluent Theory Guide equation 16.4-7.
 */
class LinearWCNSFV2PEnergyDriftFlux : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();
  LinearWCNSFV2PEnergyDriftFlux(const InputParameters & params);

  virtual Real computeElemMatrixContribution() override;

  virtual Real computeNeighborMatrixContribution() override;

  virtual Real computeElemRightHandSideContribution() override;

  virtual Real computeNeighborRightHandSideContribution() override;

  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;

  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;

  virtual void setupFaceData(const FaceInfo * face_info) override;

protected:
  /// The dimension of the simulation
  const unsigned int _dim;

  /// Dispersed phase density
  const Moose::Functor<Real> & _rho_d;

  /// Continuous phase density
  const Moose::Functor<Real> & _rho_c;

  /// Dispersed phase specific heat
  const Moose::Functor<Real> & _cp_d;

  /// Continuous phase specific heat
  const Moose::Functor<Real> & _cp_c;

  /// Dispersed phase fraction
  const Moose::Functor<Real> & _f_d;

  /// slip velocity in direction x
  const Moose::Functor<Real> & _u_slip;
  /// slip velocity in direction y
  const Moose::Functor<Real> * const _v_slip;
  /// slip velocity in direction z
  const Moose::Functor<Real> * const _w_slip;

  /// The face interpolation method for the coefficient
  const Moose::FV::InterpMethod _coeff_interp_method;

  /// Boundaries the dispersed phase may cross. The enthalpy flux carried by the relative motion
  /// is applied on every internal face and on these boundaries, and nowhere else: no phase
  /// crosses an impermeable wall, so no enthalpy may be carried through one by the relative
  /// motion, even where the temperature carries a boundary condition there.
  std::set<BoundaryID> _slip_boundaries;

  /**
   * The enthalpy flux coefficient on a face or element, beta_d beta_c / rho_m * (cp_d - cp_c).
   *
   * The phase fraction is clamped into [0, 1] before use, matching the clamping the mixture
   * property material applies.
   */
  template <typename SpaceArg>
  Real enthalpyFluxCoefficient(const SpaceArg & arg, const Moose::StateArg & state) const
  {
    const auto fd = std::clamp(_f_d(arg, state), 0.0, 1.0);
    const auto beta_d = fd * _rho_d(arg, state);
    const auto beta_c = (1.0 - fd) * _rho_c(arg, state);
    const auto rho_m = beta_d + beta_c;
    if (rho_m <= 0.0)
      return 0.0;
    return beta_d * beta_c / rho_m * (_cp_d(arg, state) - _cp_c(arg, state));
  }

private:
  /// Container for the current advected interpolation coefficients on the face, so that they are
  /// only computed once per face
  std::pair<Real, Real> _advected_interp_coeffs;

  /// The enthalpy flux on the face, cached by setupFaceData
  Real _face_flux;

  /// Multiplier that keeps the normal pointing outward on boundary faces, including boundaries
  /// which are internal to the mesh
  Real _boundary_normal_factor;

  /// The interpolation method to use for the advected temperature
  Moose::FV::InterpMethod _advected_interp_method;
};
