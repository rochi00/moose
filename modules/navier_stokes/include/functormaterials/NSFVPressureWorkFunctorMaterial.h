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
 * Computes the pressure work density carried by an advecting velocity, Dp/Dt = dp/dt + u . grad(p).
 * This is the term an energy equation written in enthalpy carries on its right hand side, and which
 * a weakly compressible formulation drops. Supplying the mixture velocity reproduces the
 * single-phase term; supplying the mixture velocity augmented by the drift contribution reproduces
 * the mixture form, in which the relative motion carries pressure work of its own.
 */
class NSFVPressureWorkFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  NSFVPressureWorkFunctorMaterial(const InputParameters & parameters);

protected:
  /// The dimension of the simulation
  const unsigned int _dim;

  /// Whether the transient part, dp/dt, is included alongside the advective part
  const bool _include_time_derivative;

  /// Pressure, the functor whose material derivative this forms
  const Moose::Functor<Real> & _pressure;

  /// Components of the velocity the pressure work is carried at
  const Moose::Functor<Real> & _u;
  const Moose::Functor<Real> * const _v;
  const Moose::Functor<Real> * const _w;
};
