//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "LinearFVElementalKernel.h"
#include "TimeIntegrator.h"

/**
 * Newton-linearized conservative enthalpy time derivative for a temperature solve.
 *
 * The current conserved quantity is temporary_density * h(T). Historical states use the accepted
 * EOS density and enthalpy. Coefficients are supplied by the configured MOOSE time integrator.
 */
class LinearFVLowMachEnthalpyTimeDerivative : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  LinearFVLowMachEnthalpyTimeDerivative(const InputParameters & params);

  Real computeMatrixContribution() override;
  Real computeRightHandSideContribution() override;

protected:
  /// Fixed end-of-step density produced by the conservative mass solve.
  const Moose::Functor<Real> & _temporary_density;

  /// Accepted EOS density, whose old state supplies rho^n.
  const Moose::Functor<Real> & _density;

  /// Exact specific enthalpy evaluated from temperature.
  const Moose::Functor<Real> & _specific_enthalpy;

  /// Current Newton derivative dh/dT.
  const Moose::Functor<Real> & _dh_dT;

  /// MOOSE time integrator supplying startup and variable-step coefficients.
  const TimeIntegrator & _time_integrator;
};
