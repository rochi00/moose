//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearWCNSFV2PDriftFluxBase.h"

/**
 * Adds the dilatation produced by the relative motion of the phases to the pressure equation of a
 * two phase mixture.
 *
 * The mass averaged velocity is not solenoidal wherever the mixture density varies, so a pressure
 * correction that enforces \f$ \nabla\cdot\mathbf{u}_m = 0 \f$ leaves out the dilatation an evolving
 * phase fraction produces. With constant phase densities and no phase change the constraint the
 * pressure equation should impose is instead
 *
 * \f[
 *   \nabla\cdot\mathbf{u}_m = -\nabla\cdot\left[\left(\alpha - c_d\right)\mathbf{u}_s\right] ,
 * \f]
 *
 * whose right hand side this object supplies, in the same form and sign convention as the
 * divergence of the predicted flux. It is an algebraic identity with no time derivative, so it holds
 * in a steady solve too, and it vanishes when the phase densities are equal since \f$ c_d = \alpha
 * \f$ then. The derivation from the volumetric flux is in the documentation page.
 */
class LinearWCNSFV2PMassDriftFlux : public LinearWCNSFV2PDriftFluxBase
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

  /// Density of the mixture
  const Moose::Functor<Real> & _rho_mixture;



  /// Cached flux for the current face, reused by the element and neighbour contributions
  Real _face_drift_flux;
};
