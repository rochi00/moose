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

class RhieChowMassFlux;

/**
 * Explicit enthalpy liquid-fraction latent heat contribution from div(rhoPhi * f_l).
 *
 * This adds the right-hand-side term
 *
 *   -L div(rhoPhi * f_l)
 *
 * to a temperature/energy equation. The supplied face flux should be a mass
 * flux per unit face area, matching the existing linear FV convention.
 */
class LinearFVLaserBeamFoamLatentHeatAdvection : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();

  LinearFVLaserBeamFoamLatentHeatAdvection(const InputParameters & params);

  Real computeElemMatrixContribution() override;
  Real computeNeighborMatrixContribution() override;
  Real computeElemRightHandSideContribution() override;
  Real computeNeighborRightHandSideContribution() override;
  Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;
  Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;

protected:
  Real computeFaceContribution();

  /// Latent heat.
  const Moose::Functor<Real> & _L;

  /// Liquid fraction.
  const Moose::Functor<Real> & _liquid_fraction;

  /// Optional Rhie-Chow mass-flux provider.
  const RhieChowMassFlux * const _mass_flux_provider;

  /// Optional directly supplied face mass-flux functor.
  const Moose::Functor<Real> * const _mass_flux_functor;

  /// Whether the directly supplied face mass flux already includes the face area.
  const bool _mass_flux_is_integrated;

  /// Cached face contribution.
  Real _face_contribution = 0.0;
};
