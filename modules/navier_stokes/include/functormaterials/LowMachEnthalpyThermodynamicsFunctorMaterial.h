//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

/**
 * Defines the exact gas/PCM enthalpy closure used by the low-Mach phase-change method.
 */
class LowMachEnthalpyThermodynamicsFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();
  LowMachEnthalpyThermodynamicsFunctorMaterial(const InputParameters & params);

protected:
  const Moose::Functor<Real> & _material_fraction;
  const Moose::Functor<Real> & _temperature;
  const Moose::Functor<Real> & _cp_gas;
  const Moose::Functor<Real> & _cp_solid;
  const Moose::Functor<Real> & _cp_liquid;
  const Moose::Functor<Real> & _rho_gas;
  const Moose::Functor<Real> & _rho_solid;
  const Moose::Functor<Real> & _rho_liquid;
  const Moose::Functor<Real> & _k_gas;
  const Moose::Functor<Real> & _k_solid;
  const Moose::Functor<Real> & _k_liquid;
  const Moose::Functor<Real> & _mu_gas;
  const Moose::Functor<Real> & _mu_solid;
  const Moose::Functor<Real> & _mu_liquid;
  const Moose::Functor<Real> & _latent_heat;
  const Moose::Functor<Real> & _T_solidus;
  const Moose::Functor<Real> & _T_liquidus;
  const Moose::Functor<Real> & _T_reference;
};
