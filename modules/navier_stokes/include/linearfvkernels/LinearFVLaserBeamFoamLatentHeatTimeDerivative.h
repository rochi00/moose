//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVElementalKernel.h"

/**
 * Explicit enthalpy liquid-fraction latent heat contribution from d(rho * f_l) / dt.
 *
 * This adds the right-hand-side term
 *
 *   -L d(rho * f_l) / dt
 *
 * to a temperature/energy equation. In LaserBeamFoam this is the time-derivative
 * part of TRHS = L * (fvc::ddt(rho, epsilon1) + fvc::div(rhoPhi, epsilon1)).
 */
class LinearFVLaserBeamFoamLatentHeatTimeDerivative : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  LinearFVLaserBeamFoamLatentHeatTimeDerivative(const InputParameters & params);

  Real computeMatrixContribution() override;
  Real computeRightHandSideContribution() override;

protected:
  /// Latent heat.
  const Moose::Functor<Real> & _L;

  /// Density multiplier.
  const Moose::Functor<Real> & _density;

  /// Optional density time derivative.
  const Moose::Functor<Real> * const _density_dot;

  /// Liquid fraction.
  const Moose::Functor<Real> & _liquid_fraction;
};
