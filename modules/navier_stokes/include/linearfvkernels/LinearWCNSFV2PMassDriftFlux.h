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

/**
 * Adds the dilatation produced by the relative motion of the phases to the pressure equation of a
 * two phase mixture.
 *
 * The mixture momentum equation is written for the mass averaged velocity \f$ \uvec_m \f$, whose
 * continuity equation is
 *
 * \f[
 *   \frac{\partial \rho_m}{\partial t} + \nabla\cdot\left(\rho_m \mathbf{u}_m\right) = 0 ,
 * \f]
 *
 * so \f$ \mathbf{u}_m \f$ is not solenoidal wherever the mixture density varies. Enforcing
 * \f$ \nabla\cdot\mathbf{u}_m = 0 \f$, as the single phase pressure correction does, therefore
 * leaves out the dilatation an evolving phase fraction produces.
 *
 * With constant phase densities and no phase change the two phase volume equations sum to
 * \f$ \nabla\cdot\mathbf{j} = 0 \f$ for the volumetric flux
 * \f$ \mathbf{j} = \alpha\mathbf{u}_d + (1-\alpha)\mathbf{u}_c \f$, which is an exact restatement
 * of the mixture continuity equation above. Using the identity
 * \f$ \mathbf{j} = \mathbf{u}_m + \left(\alpha - c_d\right)\mathbf{u}_s \f$ with
 * \f$ c_d = \alpha\rho_d/\rho_m \f$ the mass fraction of the dispersed phase, the constraint the
 * pressure equation should impose is
 *
 * \f[
 *   \nabla\cdot\mathbf{u}_m = -\nabla\cdot\left[\left(\alpha - c_d\right)\mathbf{u}_s\right] .
 * \f]
 *
 * This object adds the right hand side of that relation to the pressure equation, in the same form
 * and with the same sign convention as the divergence of the predicted flux. Being an exact
 * algebraic identity it carries no time derivative, so it is equally valid in a steady solve, where
 * \f$ \mathbf{u}_m \f$ is still not solenoidal because \f$ \rho_m \f$ varies in space.
 *
 * The contribution vanishes identically when the phase densities are equal, since \f$ c_d = \alpha
 * \f$ then, which is the check the verification test makes.
 */
class LinearWCNSFV2PMassDriftFlux : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();
  LinearWCNSFV2PMassDriftFlux(const InputParameters & params);

  virtual Real computeElemMatrixContribution() override;
  virtual Real computeNeighborMatrixContribution() override;
  virtual Real computeElemRightHandSideContribution() override;
  virtual Real computeNeighborRightHandSideContribution() override;
  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;
  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;
  virtual void setupFaceData(const FaceInfo * face_info) override;

protected:
  /// The drift contribution to the volumetric flux across the current face,
  /// \f$ (\alpha - c_d)\mathbf{u}_s\cdot\mathbf{n} \f$
  Real computeDriftVolumetricFlux();

  /// Volume fraction of the dispersed phase
  const Moose::Functor<Real> & _f_d;
  /// Density of the dispersed phase
  const Moose::Functor<Real> & _rho_d;
  /// Density of the mixture
  const Moose::Functor<Real> & _rho_mixture;

  /// Slip velocity components
  const Moose::Functor<Real> & _u_slip;
  const Moose::Functor<Real> * const _v_slip;
  const Moose::Functor<Real> * const _w_slip;

  /// Boundaries the dispersed phase may cross. The relative motion produces no dilatation across an
  /// impermeable wall, where the phases cannot separate.
  std::set<BoundaryID> _slip_boundaries;

  /// Cached flux for the current face, reused by the element and neighbour contributions
  Real _face_drift_flux;
};
