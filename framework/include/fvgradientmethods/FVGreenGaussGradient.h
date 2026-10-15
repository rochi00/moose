//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVGradientMethod.h"

#ifdef MOOSE_KOKKOS_ENABLED
namespace Moose::Kokkos
{
class GreenGaussGradient;
}
#endif

/**
 * Green-Gauss cell-centered gradient method for linear finite-volume variables.
 */
class FVGreenGaussGradient : public FVGradientMethod
{
public:
  static InputParameters validParams();
  FVGreenGaussGradient(const InputParameters & params);
  ~FVGreenGaussGradient();

private:
  void computeGradientWithoutLimiter(
      SystemBase & system,
      GradientContainer & gradient,
      const std::unordered_set<unsigned int> & variable_numbers) const override;

#ifdef MOOSE_KOKKOS_ENABLED
  void applyLimiter(SystemBase & system,
                    GradientContainer & gradient,
                    const std::unordered_set<unsigned int> & variable_numbers) const override;

  /**
   * Split variables into those whose linear FV objects are all Kokkos objects, whose gradients
   * are computed and limited with Kokkos, and the others
   * @param system The system of the variables
   * @param variable_numbers The variables
   * @returns The Kokkos variables and the others
   */
  static std::pair<std::unordered_set<unsigned int>, std::unordered_set<unsigned int>>
  splitKokkosVariables(const SystemBase & system,
                       const std::unordered_set<unsigned int> & variable_numbers);

  /// Computes the gradients of the variables whose linear FV objects are all Kokkos objects,
  /// created when first needed
  mutable std::unique_ptr<Moose::Kokkos::GreenGaussGradient> _kokkos_gradient;
#endif
};
