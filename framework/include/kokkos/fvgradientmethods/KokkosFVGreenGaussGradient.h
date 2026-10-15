//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVGradientMethod.h"

#include <memory>

namespace Moose::Kokkos
{
class GreenGaussGradient;
}

/**
 * Green-Gauss cell-centered gradient of linear finite-volume variables, computed on device (see
 * Moose::Kokkos::GreenGaussGradient). It gives the same gradients as FVGreenGaussGradient and
 * writes them into the same storage, so host and Kokkos objects can read them.
 */
class KokkosFVGreenGaussGradient : public FVGradientMethod
{
public:
  static InputParameters validParams();

  KokkosFVGreenGaussGradient(const InputParameters & params);
  ~KokkosFVGreenGaussGradient();

private:
  void computeGradientWithoutLimiter(
      SystemBase & system,
      GradientContainer & gradient,
      const std::unordered_set<unsigned int> & variable_numbers) const override;
  void applyLimiter(SystemBase & system,
                    GradientContainer & gradient,
                    const std::unordered_set<unsigned int> & variable_numbers) const override;

  /// The Kokkos computation
  const std::unique_ptr<Moose::Kokkos::GreenGaussGradient> _gradient;
};
