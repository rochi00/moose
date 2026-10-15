//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

template <typename>
class MooseLinearVariableFV;

/**
 * Test-only elemental aux kernel that outputs one component of the cell gradient of a linear FV
 * variable, computed by the variable's gradient method.
 */
class LinearFVGradientComponentTestAux : public AuxKernel
{
public:
  static InputParameters validParams();

  LinearFVGradientComponentTestAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// The linear FV variable
  MooseLinearVariableFV<Real> & _linear_variable;
  /// The gradient component
  const unsigned int _component;
};
