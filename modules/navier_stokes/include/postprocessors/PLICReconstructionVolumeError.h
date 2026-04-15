//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"
#include "MooseLinearVariableFV.h"

class LinearSystem;
class PLICReconstruction;

class PLICReconstructionVolumeError : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  PLICReconstructionVolumeError(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override {}
  virtual PostprocessorValue getValue() const override;

protected:
  const PLICReconstruction & _plic;
  MooseLinearVariableFV<Real> & _alpha_var;
  LinearSystem & _alpha_system;
  PostprocessorValue _max_alpha_error;
};
