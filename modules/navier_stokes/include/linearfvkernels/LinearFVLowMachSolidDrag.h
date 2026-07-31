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
 * Adds the low-Mach Carman-Kozeny solid penalization implicitly to momentum.
 */
class LinearFVLowMachSolidDrag : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  LinearFVLowMachSolidDrag(const InputParameters & params);

  Real computeMatrixContribution() override;
  Real computeRightHandSideContribution() override;

  bool usesCoefficient(const MooseFunctorName & coefficient_name) const
  {
    return _coefficient_name == coefficient_name;
  }

protected:
  const MooseFunctorName _coefficient_name;
  const Moose::Functor<Real> & _coefficient;
};
