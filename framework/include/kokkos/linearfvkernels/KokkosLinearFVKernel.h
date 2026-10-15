//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosLinearSystemContributionObject.h"

#include "BlockRestrictable.h"

namespace Moose::Kokkos
{

/**
 * Base class for Kokkos linear finite volume kernels that contribute to the linear system
 */
class LinearFVKernel : public LinearSystemContributionObject, public BlockRestrictable
{
public:
  static InputParameters validParams();

  LinearFVKernel(const InputParameters & parameters);
  LinearFVKernel(const LinearFVKernel & object);
};

} // namespace Moose::Kokkos
