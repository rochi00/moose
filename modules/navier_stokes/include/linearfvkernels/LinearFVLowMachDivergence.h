//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "LinearFVElementalKernel.h"

/**
 * Adds the thermodynamic low-Mach velocity-divergence constraint to the pressure equation.
 */
class LinearFVLowMachDivergence : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  LinearFVLowMachDivergence(const InputParameters & params);

  Real computeMatrixContribution() override;
  Real computeRightHandSideContribution() override;

  bool usesSource(const std::string & source_name) const { return _source_name == source_name; }

protected:
  const MooseFunctorName _source_name;
  const Moose::Functor<Real> & _source;
};
